"""Measure QC summaries of processed GSE163973 cell barcodes, without filtering."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import gc
import numpy as np
import pandas as pd
from .gse163973_reader import iter_matrices

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('archive', 'crosswalk', 'cell-metadata', 'receipt', 'output'):
        parser.add_argument('--'+option, type=Path, required=True)
    parser.add_argument('--barcode-qc-output', type=Path, help='Optional local QC metadata; keep outside source files and version control')
    args = parser.parse_args()
    inputs=[args.archive,args.crosswalk,args.cell_metadata]
    outputs = [args.output] + ([args.barcode_qc_output] if args.barcode_qc_output else [])
    if any(p.exists() for p in outputs) or len({p.resolve() for p in outputs}) != len(outputs):
        raise ValueError('Outputs must be distinct new files; no source or existing file is overwritten')
    before = {name: sha(path) for name,path in zip(('archive','crosswalk','cell_metadata'), inputs)}
    receipt = json.loads(args.receipt.read_text())
    expected = {'archive': receipt['source_files']['GSE163973_RAW.tar']['sha256'],
                'crosswalk': receipt['output_sha256']['cell_barcode_crosswalk.parquet'],
                'cell_metadata': receipt['output_sha256']['source_cell_metadata.parquet']}
    if before != expected:
        raise ValueError('A processed source does not match the supplied recovery receipt')
    cross=pd.read_parquet(inputs[1])
    all_cells=[]
    for sample, genes, barcodes, matrix in iter_matrices(args.archive, strict_members=True, thread_limit=2):
        print('Measuring',sample,flush=True)
        sub=cross.loc[cross.source_donor_key.eq(sample)].set_index('matrix_barcode',verify_integrity=True).loc[barcodes]
        counts=np.asarray(matrix.sum(axis=0)).ravel()
        if not np.array_equal(counts,sub.matrix_sum_UMI):
            raise ValueError('Matrix totals differ from barcode crosswalk')
        detected=np.diff(matrix.indptr)
        mito=genes.iloc[:, 1].str.upper().str.startswith('MT-').to_numpy()
        mt=np.asarray(matrix[mito].sum(axis=0)).ravel()
        frac=np.divide(mt*100,counts,out=np.full(len(counts),np.nan),where=counts>0)
        d=pd.DataFrame({'sample':sample,'barcode':barcodes,'annotated':sub.source_metadata_row.notna().to_numpy(),
            'cell_type':sub.cellType.fillna('unannotated').to_numpy(),'matrix_counts':counts,
            'matrix_detected_features':detected,'matrix_mito_pct':frac,
            'author_counts':sub.nCount_RNA.to_numpy(),'author_detected_features':sub.nFeature_RNA.to_numpy(),
            'author_mito_pct':sub['percent.mt'].to_numpy()})
        all_cells.append(d)
        del matrix
        gc.collect()
    cells=pd.concat(all_cells,ignore_index=True)
    if len(cells) != 53148 or cells.annotated.sum() != 40655 or cells.duplicated(['sample','barcode']).any():
        raise ValueError('Unexpected barcode count, annotations, or duplicate barcodes')
    # MAD flags are descriptive review markers; no cells are excluded.
    cells['flag_low_counts']=False; cells['flag_low_features']=False; cells['flag_high_mito']=False
    thresholds=[]
    for sample,idx in cells.loc[cells.annotated].groupby('sample').groups.items():
        for col,flag,transform,direction in [('matrix_counts','flag_low_counts','log1p','low'),('matrix_detected_features','flag_low_features','log1p','low'),('matrix_mito_pct','flag_high_mito','none','high')]:
            x=cells.loc[idx,col].to_numpy(float)
            v=np.log1p(x) if transform=='log1p' else x
            med=np.nanmedian(v); mad=np.nanmedian(abs(v-med)); threshold=med+(-3 if direction=='low' else 3)*mad
            flagvalues=v<threshold if direction=='low' else v>threshold
            cells.loc[idx,flag]=flagvalues
            thresholds.append({'sample':sample,'metric':col,'transform':transform,'median':float(med),'mad':float(mad),'threshold':float(threshold),'direction':direction,'n_flagged':int(flagvalues.sum())})
    cells['flag_any']=cells[['flag_low_counts','flag_low_features','flag_high_mito']].any(axis=1)
    def q(x):
        a=x.to_numpy(dtype=float); return {f'q{int(k*100):02d}':float(np.nanquantile(a,k)) for k in [.01,.05,.5,.95,.99]}
    summaries=[]
    for (sample,annotated),d in cells.groupby(['sample','annotated']):
        row={'sample':sample,'annotated':bool(annotated),'barcodes':len(d),'counts':q(d.matrix_counts),'detected':q(d.matrix_detected_features),'mito_pct':q(d.matrix_mito_pct),'counts_over_6000':int((d.matrix_counts>6000).sum()),'mito_over_15pct':int((d.matrix_mito_pct>15).sum()),'review_flag_any':int(d.flag_any.sum())}
        if annotated:
            row['author_feature_exact_fraction']=float((d.matrix_detected_features==d.author_detected_features).mean())
            row['author_feature_absdiff_q99']=float((d.matrix_detected_features-d.author_detected_features).abs().quantile(.99))
        summaries.append(row)
    types=[]
    for (sample,typ),d in cells.loc[cells.annotated].groupby(['sample','cell_type']):
        types.append({'sample':sample,'cell_type':typ,'cells':len(d),'count_median':float(d.matrix_counts.median()),'features_median':float(d.matrix_detected_features.median()),'mito_median':float(d.matrix_mito_pct.median()),'review_flag_any':int(d.flag_any.sum())})
    # QC metadata is separate from, and does not rewrite, expression tables.
    if args.barcode_qc_output:
        args.barcode_qc_output.parent.mkdir(parents=True, exist_ok=True)
        cells.to_parquet(args.barcode_qc_output,index=False)
    if before != {name:sha(path) for name,path in zip(('archive','crosswalk','cell_metadata'), inputs)}:
        raise ValueError('An input changed during QC')
    result={'samples':summaries,'cell_types':types,'descriptive_thresholds':thresholds,
       'total_barcodes':len(cells),'annotated_barcodes':int(cells.annotated.sum()),
       'annotated_counts_over_6000':int((cells.annotated & cells.matrix_counts.gt(6000)).sum()),
       'annotated_review_flag_any':int((cells.annotated & cells.flag_any).sum()),
       'source_unchanged':True,'input_sha256':before,'script_sha256':sha(Path(__file__)),
       'reader_sha256':sha(Path(iter_matrices.__code__.co_filename)),
       'notes':['Three unscaled median absolute deviations are descriptive only, per sample on annotated cells; do not use to remove cells without review.',
                'Count and mitochondrial fields use complete processed matrices; exclude from model input and do not use as hidden-gene-derived normalization.',
                'No unannotated barcode is declared a valid cell; no doublet scoring or new expression matrix produced.']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['total_barcodes','annotated_barcodes','annotated_counts_over_6000','annotated_review_flag_any','source_unchanged']},indent=2),flush=True)

if __name__=='__main__': main()
