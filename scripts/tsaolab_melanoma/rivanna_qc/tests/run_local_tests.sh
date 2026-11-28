#!/bin/bash
# Local tests for the rivanna_qc quantification scripts, using synthetic fixtures and stub executables
# (salmon, apptainer/fastp, cutadapt, fastqc, bbversion.sh, and an Lmod-like `module` function that reproduces the
# observed bzip2/1.0.6 -> 1.0.8 swap when Apptainer is loaded). No network, Slurm, real data or real tools.
# These tests exercise control flow, fail-closed behaviour and provenance only; they say nothing about the
# numerical behaviour of the real fastp/Salmon binaries on Rivanna.
#
# Usage: bash scripts/tsaolab_melanoma/rivanna_qc/tests/run_local_tests.sh [python3]
#        bash .../run_local_tests.sh --isolation-probe [python3]   (print the isolated environment and exit)
#
# Self-isolation: the suite re-executes itself under `env -i ... bash --noprofile --norc` with a minimal PATH, so
# nothing from the caller's shell reaches the stubs. On an HPC login node this matters: Lmod exports a `module`
# shell function (BASH_FUNC_module%%), LMOD_*/MODULEPATH/_ModuleTable state, and BASH_ENV pointing at Lmod's
# init script, which every non-interactive child `bash` re-sources -- silently replacing the stub `module` with
# the real one (the cause of 13 failures observed on Rivanna). Production scripts are unaffected: they are meant
# to use the real Lmod and are not changed by this. Set RQC_TESTS_NO_ISOLATE=1 only to debug the harness itself.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
QC="$(cd "$HERE/.." && pwd)"
PROBE=0; [ "${1:-}" = --isolation-probe ] && { PROBE=1; shift; }
PY=${1:-python3}
if [ "${RQC_TESTS_ISOLATED:-0}" != 1 ] && [ "${RQC_TESTS_NO_ISOLATE:-0}" != 1 ]; then
  PYABS=$(command -v "$PY" 2>/dev/null) || { echo "ERROR: python interpreter '$PY' not found" >&2; exit 2; }
  case "$PYABS" in /*) ;; *) PYABS="$(cd "$(dirname "$PYABS")" && pwd)/$(basename "$PYABS")" ;; esac
  BASH_ABS=$(command -v bash)
  args=(); [ "$PROBE" = 1 ] && args+=(--isolation-probe)
  exec env -i HOME="${HOME:-/tmp}" TMPDIR="${TMPDIR:-/tmp}" LANG=C LC_ALL=C TERM=dumb \
    PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin" RQC_TESTS_ISOLATED=1 KEEP_TEST_DIR="${KEEP_TEST_DIR:-}" \
    "$BASH_ABS" --noprofile --norc "$0" ${args[@]+"${args[@]}"} "$PYABS"
fi
# Second line of defence (also covers RQC_TESTS_NO_ISOLATE=1): drop inherited functions and module state.
unset BASH_ENV ENV
for _fn in $(compgen -A function); do unset -f "$_fn"; done
for _v in $(compgen -e); do case "$_v" in LMOD*|__LMOD*|_ModuleTable*|MODULEPATH*|MODULESHOME|LOADEDMODULES|_LMFILES_|__Init_Default_Modules|ML_*) unset "$_v" ;; esac; done
unset _fn _v
if [ "$PROBE" = 1 ]; then
  echo "isolated=${RQC_TESTS_ISOLATED:-0}"; echo "BASH_ENV=${BASH_ENV:-<unset>}"
  echo "module=$(type -t module 2>/dev/null || echo none)"; echo "child_module=$(bash -c 'type -t module' 2>/dev/null || echo none)"
  echo "lmod_vars=$(compgen -e | grep -cE '^(LMOD|__LMOD|_ModuleTable|MODULEPATH|LOADEDMODULES)' || true)"
  echo "path=$PATH"; exit 0
fi
T=$(mktemp -d "${TMPDIR:-/tmp}/rqc_tests.XXXXXX")
trap "chmod -R u+w \"$T\" 2>/dev/null; [ -n \"${KEEP_TEST_DIR:-}\" ] || rm -rf \"$T\"" EXIT
pass=0; fail=0
ok()   { local n=$1; shift; if "$@" >"$T/last.log" 2>&1; then pass=$((pass+1)); echo "PASS  $n"; else fail=$((fail+1)); echo "FAIL  $n (expected success)"; sed 's/^/      /' "$T/last.log" | tail -15; fi; }
bad()  {  # bad NAME EXPECTED_MESSAGE_REGEX CMD...: must fail AND fail for the expected reason
  local n=$1 why=$2; shift 2
  if "$@" >"$T/last.log" 2>&1; then fail=$((fail+1)); echo "FAIL  $n (expected failure)"; sed 's/^/      /' "$T/last.log" | tail -15
  elif ! grep -Eq -- "$why" "$T/last.log"; then fail=$((fail+1)); echo "FAIL  $n (failed, but not with /$why/)"; sed 's/^/      /' "$T/last.log" | tail -8
  else pass=$((pass+1)); echo "PASS  $n  [$(grep -Eo -- "$why" "$T/last.log" | head -1 | cut -c1-90)]"; fi; }
# has/hasnt FAIL when the file is missing or unreadable: a negative assertion must never pass on absent evidence.
has()  { local n=$1 pat=$2 f=$3; if [ ! -f "$f" ] || [ ! -r "$f" ]; then fail=$((fail+1)); echo "FAIL  $n (file missing: $f)"
         elif grep -Eq -- "$pat" "$f"; then pass=$((pass+1)); echo "PASS  $n"; else fail=$((fail+1)); echo "FAIL  $n (/$pat/ not in $f)"; fi; }
hasnt(){ local n=$1 pat=$2 f=$3; if [ ! -f "$f" ] || [ ! -r "$f" ]; then fail=$((fail+1)); echo "FAIL  $n (file missing: $f)"
         elif grep -Eq -- "$pat" "$f"; then fail=$((fail+1)); echo "FAIL  $n (/$pat/ found in $f)"; else pass=$((pass+1)); echo "PASS  $n"; fi; }
selftest() ( pass=0; fail=0; "$@" >/dev/null 2>&1; [ "$fail" -eq 1 ] && [ "$pass" -eq 0 ] )  # helper must record a failure
cond() { local n=$1; shift; if "$@"; then pass=$((pass+1)); echo "PASS  $n"; else fail=$((fail+1)); echo "FAIL  $n"; fi; }
flat() { if [ -f "$1" ]; then tr -d '\n' < "$1" > "$1.flat"; echo "$1.flat"; else echo "$1.missing"; fi; }
h256() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
hmd5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
TOOLS="$QC/rivanna_qc_tools.py"

echo "== harness self-tests"
cond "has fails on a missing file"              selftest has x . "$T/does-not-exist"
cond "hasnt fails on a missing file"            selftest hasnt x . "$T/does-not-exist"
cond "bad fails when the reason does not match" selftest bad x "no-such-message" false
cond "bad fails when the command succeeds"      selftest bad x . true
cond "suite runs isolated (env -i)"             test "${RQC_TESTS_ISOLATED:-0}" = 1 -o "${RQC_TESTS_NO_ISOLATE:-0}" = 1
cond "no BASH_ENV inherited"                    test -z "${BASH_ENV:-}"
cond "no Lmod state inherited"                  test -z "$(compgen -e | grep -E '^(LMOD|__LMOD|_ModuleTable|MODULEPATH|LOADEDMODULES)')"

echo "== static checks"
for f in "$QC"/*.sh "$QC"/*.sbatch "$QC"/tests/*.sh; do ok "bash -n $(basename "$f")" bash -n "$f"; done
ok "python parses rivanna_qc_tools.py" "$PY" -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$TOOLS"
if command -v shellcheck >/dev/null; then ok "shellcheck" shellcheck -S warning "$QC"/*.sh "$QC"/*.sbatch; else echo "SKIP  shellcheck (not installed)"; fi

# ------------------------------------------------------------------------------------------------ stubs
STUB="$T/bin"; mkdir -p "$STUB"
cat > "$STUB/salmon" <<'EOF'
#!/bin/bash
set -e
case "$1" in
  --version) echo "salmon ${STUB_SALMON_VERSION:-1.10.3}" ;;
  index) shift; while [ $# -gt 0 ]; do case "$1" in -t) fa=$2; shift;; -i) idx=$2; shift;; -k) k=$2; shift;; esac; shift; done
    gzip -t "$fa"; mkdir "$idx"
    cat > "$idx/info.json" <<J
{"SeqHash": "${STUB_SEQ_HASH}", "NameHash": "${STUB_NAME_HASH}", "SeqHash512": "${STUB_SEQ_HASH512}",
 "NameHash512": "${STUB_NAME_HASH512}", "DecoySeqHash": "${STUB_DECOY_HASH}", "DecoyNameHash": "${STUB_DECOY_HASH}", "k": $k}
J
    echo '{"indexVersion": 2}' > "$idx/versionInfo.json" ;;
  quant) shift; while [ $# -gt 0 ]; do case "$1" in -i) idx=$2; shift;; -1) r1=$2; shift;; -2) r2=$2; shift;; -o) out=$2; shift;; esac; shift; done
    gzip -t "$r1"; gzip -t "$r2"; mkdir -p "$out/aux_info"
    n1=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
    if [ "${STUB_TRUNCATE:-0}" = 1 ]; then head -100 "$STUB_QUANT_SF" > "$out/quant.sf"; else cp "$STUB_QUANT_SF" "$out/quant.sf"; fi
    rows=$(( $(wc -l < "$STUB_QUANT_SF") - 1 ))
    "$STUB_PY" - "$idx/info.json" "$out/aux_info/meta_info.json" "${STUB_PCT:-35.02}" "$rows" "$n1" <<'P'
import json,sys
i=json.load(open(sys.argv[1])); n=int(sys.argv[5])
m={"salmon_version":"1.10.3","percent_mapped":float(sys.argv[3]),"num_valid_targets":int(sys.argv[4]),
   "num_processed":n,"num_mapped":int(n*float(sys.argv[3])/100),
   "index_seq_hash":i["SeqHash"],"index_name_hash":i["NameHash"],"index_seq_hash512":i["SeqHash512"],
   "index_name_hash512":i["NameHash512"],"index_decoy_seq_hash":i["DecoySeqHash"],"index_decoy_name_hash":i["DecoyNameHash"]}
json.dump(m,open(sys.argv[2],"w"))
P
    echo "{\"expected_format\": \"${STUB_LIBTYPE:-ISR}\", \"compatible_fragment_ratio\": 1.0}" > "$out/lib_format_counts.json" ;;
  *) echo "stub salmon: unsupported $*" >&2; exit 2 ;;
esac
EOF
cat > "$STUB/apptainer" <<'EOF'
#!/bin/bash
set -e
case "$1" in
  --version) echo "apptainer version 1.4.5" ;;
  pull) printf '%s' "${STUB_SIF_CONTENT:-stub-fastp-sif-0.23.4}" > "$2" ;;
  exec) shift 2; [ "$1" = fastp ] || exit 2; shift
    if [ "${1:-}" = --version ]; then echo "fastp ${STUB_FASTP_VERSION:-0.23.4}" >&2; exit 0; fi
    while [ $# -gt 0 ]; do case "$1" in -i) i1=$2; shift;; -I) i2=$2; shift;; -o) o1=$2; shift;; -O) o2=$2; shift;; -j) j=$2; shift;; -h) h=$2; shift;; esac; shift; done
    cp "$i1" "$o1"; cp "$i2" "$o2"
    n=$(( ($(gzip -dc "$i1" | wc -l) + $(gzip -dc "$i2" | wc -l)) / 4 ))
    echo "{\"summary\":{\"fastp_version\":\"0.23.4\",\"before_filtering\":{\"total_reads\":$n},\"after_filtering\":{\"total_reads\":$n}}}" > "$j"
    echo '<html></html>' > "$h" ;;
  *) exit 2 ;;
esac
EOF
printf '#!/bin/bash\necho 5.2\n' > "$STUB/cutadapt"
printf '#!/bin/bash\necho "FastQC v0.12.1"\n' > "$STUB/fastqc"
printf '#!/bin/bash\necho 39.01\n' > "$STUB/bbversion.sh"
chmod +x "$STUB"/*
# Lmod-like stub: salmon pulls bzip2/1.0.6; loading apptainer swaps it to bzip2/1.0.8 (as observed on Rivanna).
# Every call is appended to $RQC_STUB_MODULE_LOG so tests can prove the stub (not a real Lmod) was invoked.
module() {
  : RQC_STUB_MODULE
  local m
  echo "$*" >> "${RQC_STUB_MODULE_LOG:-/dev/null}"
  case "$1" in
    purge) _STUB_MODS="" ;;
    load) shift; for m in "$@"; do
        case "$m" in
          salmon/*) _STUB_MODS="$_STUB_MODS bzip2/1.0.6 $m" ;;
          apptainer*) _STUB_MODS="${_STUB_MODS//bzip2\/1.0.6/} bzip2/1.0.8 apptainer/1.4.5" ;;
          nonexistent*) echo "Lmod has detected the following error: $m" >&2; return 1 ;;
          *) _STUB_MODS="$_STUB_MODS $m" ;;
        esac; done ;;
    -t) for m in ${_STUB_MODS:-}; do echo "$m" >&2; done ;;
  esac
}
export -f module
export RQC_STUB_MODULE_LOG="$T/module_calls.log"; : > "$RQC_STUB_MODULE_LOG"
stub_module_reaches_child() {  # ${@+} for bash 3.2 + set -u; output captured first (grep -q + pipefail races SIGPIPE)
  local t; t=$(env ${@+"$@"} bash -c 'type module 2>/dev/null' || true); [[ "$t" == *RQC_STUB_MODULE* ]]; }
# A login shell's Lmod, reproduced: BASH_ENV re-defines `module` in every child bash and puts "real" tools first.
FOREIGN="$T/foreign_lmod"; mkdir -p "$FOREIGN/bin"
printf '#!/bin/bash\n[ "$1" = --version ] && { echo "apptainer version 1.4.5"; exit 0; }\necho "FATAL:   could not open image $2: image format not recognized" >&2; exit 255\n' > "$FOREIGN/bin/apptainer"
printf '#!/bin/bash\n[ "$1" = --version ] && { echo "salmon 1.10.3"; exit 0; }\necho "real salmon: cannot read stub index" >&2; exit 1\n' > "$FOREIGN/bin/salmon"
chmod +x "$FOREIGN/bin/"*
cat > "$FOREIGN/init.sh" <<EOF
module() { echo "FOREIGN \$*" >> "$FOREIGN/calls.log"; case "\$1" in load) export PATH="$FOREIGN/bin:\$PATH";; -t) echo gcc/14.2.0 >&2;; esac; }
EOF

# --------------------------------------------------------------------------------------------- fixtures
F="$T/fixtures"; mkdir -p "$F/gencode" "$F/raw"
printf '>ENST00000000001.1|ENSG00000000001.1|-|-|TX1-201|TX1|8|protein_coding|\nACGTACGT\n' | gzip -n > "$F/gencode/gencode.v46.transcripts.fa.gz"
GMD5=$(hmd5 "$F/gencode/gencode.v46.transcripts.fa.gz")
printf '%s  gencode.v46.transcripts.fa.gz\n0123456789abcdef0123456789abcdef  other.gtf.gz\n' "$GMD5" > "$F/gencode/MD5SUMS"
SPECS="01 02 03 37 43"
for s in $SPECS; do for r in 1 2; do
  printf '@r1/%s\nACGTACGTAGATCGGAAG\n+\nIIIIIIIIIIIIIIIIII\n@r2/%s\nGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIII\n@r3/%s\nACGTTTTTACGTTTTTAC\n+\nIIIIIIIIIIIIIIIIII\n' "$r" "$r" "$r" | gzip -n > "$F/raw/20088XR-11-${s}_S0_L001_R${r}_001.fastq.gz"
done; done
{ echo "sample_id,core_id,project,species,tissue,specimen,disease_label,lane,index_s,r1,r2,md5_r1,md5_r2"
  for s in $SPECS; do echo "GM$s,20088XR-11-$s,20088-11,human,melanoma,FFPE,melanoma,L001,S0,20088XR-11-${s}_S0_L001_R1_001.fastq.gz,20088XR-11-${s}_S0_L001_R2_001.fastq.gz,$(hmd5 "$F/raw/20088XR-11-${s}_S0_L001_R1_001.fastq.gz"),$(hmd5 "$F/raw/20088XR-11-${s}_S0_L001_R2_001.fastq.gz")"; done; } > "$F/samples.csv"
# validated pairs: 01 (core-id label), 02 (short label), 37 (gzip-valid but blocked); 03 absent; 43 absent (corrupt)
{ printf 'sample\tr1\tr2\n'
  printf '20088XR-11-01\t%s\t%s\n' "$F/raw/20088XR-11-01_S0_L001_R1_001.fastq.gz" "$F/raw/20088XR-11-01_S0_L001_R2_001.fastq.gz"
  printf '11-02\t%s\t%s\n' "$F/raw/20088XR-11-02_S0_L001_R1_001.fastq.gz" "$F/raw/20088XR-11-02_S0_L001_R2_001.fastq.gz"
  printf '20088XR-11-37\t%s\t%s\n' "$F/raw/20088XR-11-37_S0_L001_R1_001.fastq.gz" "$F/raw/20088XR-11-37_S0_L001_R2_001.fastq.gz"; } > "$F/validated_pairs.tsv"
"$PY" - "$F" <<'P'
import random,sys,gzip
random.seed(1); F=sys.argv[1]
rows=[(f"ENST{i:011d}.1",1000,900.0,random.expovariate(0.05) if i%4 else 0.0) for i in range(2000)]
s=sum(r[3] for r in rows)
with open(f"{F}/reference_quant.sf","w") as fh:
    fh.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
    for n,l,e,t in rows: fh.write(f"{n}\t{l}\t{e}\t{t/s*1e6:.6f}\t{t:.3f}\n")
with gzip.open(f"{F}/genemap.fa.gz","wt") as o:
    for n,*_ in rows: o.write(f">{n}|ENSG{int(n[4:15])//3:011d}.1|-|-|x|x|1000|protein_coding|\nACGT\n")
P
chmod -R a-w "$F/raw"
RAW_BEFORE=$(for f in "$F"/raw/*; do h256 "$f"; done | tr '\n' ' ')

# ---------------------------------------------------------------------------------------- environment
export PATH="$STUB:$PATH" RQC_SKIP_MODULE_LOAD=1 PYTHON_BIN="$PY" RQC_DIR="$QC"
export EXPECT_SEQ_HASH=9edf59907918a715a58542b94221f7ad7d88c2edf40104a2ce1834a37c9df644
export EXPECT_NAME_HASH=11bd2542a03abacba2f4ba1e5f2a8b8a8402fd88fc3e55050569893fa17b4b94
export EXPECT_SEQ_HASH512=7eb64cfebb1fa1ea843da8e5a570842a873c09f177bd57815fd01390d5a1ae56ffef711b9100383eb57cb786473d5c45baed292699bc91cadeaeb498f2446a19
export EXPECT_NAME_HASH512=b22b5cb3f38f0a919fc7a2a834f044c92b90ffc19ade1e4abc40ecf520b2b58ffd7b4f57a73028fe0c4a83b01bef6c591d113aa3f662d7ad620f72dca0d5ac3d
export EXPECT_DECOY_HASH=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
export STUB_SEQ_HASH=$EXPECT_SEQ_HASH STUB_NAME_HASH=$EXPECT_NAME_HASH STUB_SEQ_HASH512=$EXPECT_SEQ_HASH512 \
       STUB_NAME_HASH512=$EXPECT_NAME_HASH512 STUB_DECOY_HASH=$EXPECT_DECOY_HASH STUB_QUANT_SF="$F/reference_quant.sf" STUB_PY="$PY"
W="$T/work"; mkdir -p "$W"
SIF_SHA=$(printf '%s' stub-fastp-sif-0.23.4 | { command -v sha256sum >/dev/null && sha256sum || shasum -a 256; } | cut -d' ' -f1)
export PROVENANCE_DIR="$W/provenance" REF_DIR="$W/ref" GENCODE_BASE_URL="file://$F/gencode" \
       GENCODE_TRANSCRIPTS=gencode.v46.transcripts.fa.gz GENCODE_TRANSCRIPTS_MD5=$GMD5 \
       INDEX_DIR="$W/ref/salmon_index_k31" SALMON_KMER=31 \
       FASTP_IMAGE_URI=docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5 FASTP_SIF="$W/containers/fastp.sif" \
       FASTP_EXPECTED_VERSION=0.23.4 FASTP_SIF_SHA256=$SIF_SHA \
       SAMPLE_SHEET="$F/samples.csv" VALIDATED_PAIRS="$F/validated_pairs.tsv" RAW_DIR="$F/raw" QUANT_ROOT="$W/quant" \
       EXPECT_NUM_TARGETS=2000 EXPECT_LIBTYPE=ISR EXCLUDED_SPECIMENS="" LIFT_EXCLUSIONS="" \
       SALMON_MODULES="gcc/14.2.0 salmon/1.10.3" APPTAINER_MODULES="apptainer"
COMMITTED_META="$T/committed_GM01_meta_info.json"
git -C "$QC" show HEAD:data/processed/tsaolab_melanoma/salmon_quant/GM01/aux_info/meta_info.json > "$COMMITTED_META" 2>/dev/null || rm -f "$COMMITTED_META"
Q="$QC/quantify_specimen.sbatch"

echo "== verify_environment.sh"
ok  "environment passes with stubs"            bash "$QC/verify_environment.sh" "$W/env1.tsv"
has "environment report lists bbmap PASS"       $'bbmap\t.*PASS' "$W/env1.tsv"
bad "environment fails on wrong cutadapt"       "FAIL  cutadapt" env EXPECT_CUTADAPT=9.9 bash "$QC/verify_environment.sh" "$W/env2.tsv"
bad "environment refuses report collision"      "refusing to overwrite" bash "$QC/verify_environment.sh" "$W/env1.tsv"
echo "== inherited-environment isolation (module functions, BASH_ENV, Lmod state)"
cond "stub module reaches child bash"           stub_module_reaches_child
cond "foreign BASH_ENV would replace it (guard detects)" bash -c 't=$(env BASH_ENV="$1" bash -c "type module"); [[ "$t" != *RQC_STUB_MODULE* && "$t" == *FOREIGN* ]]' _ "$FOREIGN/init.sh"
FMOD='() {  echo FOREIGN-EXPORTED "$@" >&2; return 1
}'
env -i PATH=/usr/bin:/bin "BASH_FUNC_module%%=$FMOD" BASH_ENV="$FOREIGN/init.sh" LMOD_CMD=/x/lmod MODULEPATH=/x \
  bash "$HERE/run_local_tests.sh" --isolation-probe "$PY" > "$W/probe.txt" 2>&1
has "contaminated caller: suite re-execs isolated" '^isolated=1$' "$W/probe.txt"
has "  ... inherited BASH_ENV dropped"            '^BASH_ENV=<unset>$' "$W/probe.txt"
has "  ... inherited exported module() dropped"   '^module=none$' "$W/probe.txt"
has "  ... child bash sees no module()"           '^child_module=none$' "$W/probe.txt"
has "  ... Lmod variables dropped"                '^lmod_vars=0$' "$W/probe.txt"
ok  "environment with Lmod stub"               env RQC_SKIP_MODULE_LOAD=0 bash "$QC/verify_environment.sh" "$W/env3.tsv"
has "  ... stub module was actually invoked"     '^load cutadapt/5\.2$' "$RQC_STUB_MODULE_LOG"
has "  ... salmon modules loaded via stub"       '^load gcc/14\.2\.0 salmon/1\.10\.3$' "$RQC_STUB_MODULE_LOG"
has "co-load swap bzip2 1.0.6->1.0.8 reported"  'coload_swaps.*bzip2/1\.0\.6->bzip2/1\.0\.8' "$W/env3.tsv"
bad "unloadable module set fails"               "FAIL  salmon" env RQC_SKIP_MODULE_LOAD=0 SALMON_MODULES=nonexistent/1 bash "$QC/verify_environment.sh" "$W/env4.tsv"
has "  ... failure came from the stub module"    '^load nonexistent/1$' "$RQC_STUB_MODULE_LOG"

echo "== fetch_gencode_v46.sh"
ok  "fetch downloads and verifies"             bash "$QC/fetch_gencode_v46.sh"
ok  "fetch rerun verifies existing (no overwrite)" bash "$QC/fetch_gencode_v46.sh"
bad "fetch fails when pinned MD5 disagrees"     "disagrees with pinned" env GENCODE_TRANSCRIPTS_MD5=00000000000000000000000000000000 bash "$QC/fetch_gencode_v46.sh"
bad "fetch rejects malformed MD5"               "not a 32-hex MD5" env GENCODE_TRANSCRIPTS_MD5=xyz bash "$QC/fetch_gencode_v46.sh"
mkdir -p "$W/ref_bad"; cp "$F/gencode/MD5SUMS" "$W/ref_bad/MD5SUMS.release_46"; printf 'corrupt' > "$W/ref_bad/gencode.v46.transcripts.fa.gz"
bad "fetch fails on corrupt existing file"      "MD5 mismatch for" env REF_DIR="$W/ref_bad" bash "$QC/fetch_gencode_v46.sh"
has "corrupt existing file left untouched"      corrupt "$W/ref_bad/gencode.v46.transcripts.fa.gz"
mkdir -p "$W/ref_link"; cp "$F/gencode/MD5SUMS" "$W/ref_link/MD5SUMS.release_46"; ln -s "$W/nowhere" "$W/ref_link/gencode.v46.transcripts.fa.gz"
bad "fetch will not adopt over a dangling symlink" "cannot adopt" env REF_DIR="$W/ref_link" bash "$QC/fetch_gencode_v46.sh"
cond "dangling symlink target never written"    test ! -e "$W/nowhere"
cond "no partial download left behind"          test ! -e "$REF_DIR/gencode.v46.transcripts.fa.gz.partial"

echo "== provision_fastp_container.sh"
bad "provision rejects tag not pinning version" "does not pin fastp" env FASTP_IMAGE_URI=docker://quay.io/biocontainers/fastp:latest bash "$QC/provision_fastp_container.sh"
bad "provision rejects digest not in URI"       "must pin @sha256" env FASTP_IMAGE_DIGEST=sha256:$(printf '%064d' 1) bash "$QC/provision_fastp_container.sh"
bad "provision fails closed on new SIF bytes"   "SIF SHA-256 mismatch" env STUB_SIF_CONTENT=different bash "$QC/provision_fastp_container.sh"
has "  ... and explains independent pulls differ" "Independently pulled SIFs can differ" "$T/last.log"
cond "mismatched pull not adopted"              test ! -e "$FASTP_SIF"
rm -f "$FASTP_SIF.partial"
ok  "provision pulls and verifies"             bash "$QC/provision_fastp_container.sh"
ok  "provision rerun verifies existing"        bash "$QC/provision_fastp_container.sh"
bad "provision fails on wrong fastp version"    "unexpected fastp version" env STUB_FASTP_VERSION=0.24.0 bash "$QC/provision_fastp_container.sh"
ok  "new artifact adopted only with explicit opt-in" env STUB_SIF_CONTENT=rebuilt FASTP_ACCEPT_NEW_SIF=1 FASTP_SIF="$W/containers/fastp_new.sif" bash "$QC/provision_fastp_container.sh"
has "  ... and asks to update pinned hash"      "ACTION: set FASTP_SIF_SHA256=" "$T/last.log"
has "  ... provenance records new-artifact mode" "pulled_new_artifact_accepted" "$(ls -t "$PROVENANCE_DIR"/provision_fastp_*.json | head -1)"
bad "opt-in still requires fastp version"       "not adopted" env STUB_SIF_CONTENT=rebuilt2 STUB_FASTP_VERSION=0.20.0 FASTP_ACCEPT_NEW_SIF=1 FASTP_SIF="$W/containers/fastp_new2.sif" bash "$QC/provision_fastp_container.sh"

echo "== build_salmon_index.sbatch"
bad "index build refuses fingerprint mismatch"  "did not all MATCH" env STUB_SEQ_HASH=deadbeef INDEX_DIR="$W/ref/idx_bad" bash "$QC/build_salmon_index.sbatch"
cond "mismatched index not adopted"             test ! -e "$W/ref/idx_bad"
cond "lock released after failure"              test ! -e "$W/ref/idx_bad.lock"
mkdir "$W/ref/idx_locked.lock"
bad "index build blocked by concurrent lock"    "lock held or stale" env INDEX_DIR="$W/ref/idx_locked" bash "$QC/build_salmon_index.sbatch"
ok  "index build verifies six fingerprints + k" bash "$QC/build_salmon_index.sbatch"
bad "index build refuses existing index"        "refusing to overwrite" bash "$QC/build_salmon_index.sbatch"
bad "index build fails on missing transcriptome" "missing, empty or unreadable" env GENCODE_TRANSCRIPTS=missing.fa.gz INDEX_DIR="$W/ref/idx2" bash "$QC/build_salmon_index.sbatch"
BP=$(flat "$(ls "$PROVENANCE_DIR"/build_salmon_index_*.json | head -1)")
has "index provenance: k is candidate"          '"kmer": \{[^}]*"status": "candidate"' "$BP"
has "index provenance: decoys verified-collaborator" '"decoys": \{[^}]*"status": "verified-collaborator"' "$BP"

echo "== verify-fingerprint"
ok  "index dir: all six MATCH"                 bash "$QC/verify_reference_fingerprint.sh" "$INDEX_DIR" "$W/fp1.json"
has "report records collaborator k unverified"  "UNVERIFIED" "$W/fp1.json"
if [ -s "$COMMITTED_META" ]; then
  ok  "committed GM01 meta_info: all six MATCH (incl. SHA-512)" bash "$QC/verify_reference_fingerprint.sh" "$COMMITTED_META" "$W/fp_committed.json"
  has "  ... seq_hash512 and name_hash512 MATCH" '"seq_hash512",[[:space:]]*"status": "MATCH".*"name_hash512",[[:space:]]*"status": "MATCH"' "$(flat "$W/fp_committed.json")"
else echo "SKIP  committed GM01 meta_info (blob not available)"; fi
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); d['NameHash512']='00'; json.dump(d,open(sys.argv[2],'w'))" "$INDEX_DIR/info.json" "$W/info_mm.json"
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); del d['SeqHash512']; json.dump(d,open(sys.argv[2],'w'))" "$INDEX_DIR/info.json" "$W/info_na.json"
bad "NameHash512 mismatch detected"             "name_hash512 +MISMATCH" bash "$QC/verify_reference_fingerprint.sh" "$W/info_mm.json" "$W/fp2.json"
bad "missing SeqHash512 -> UNAVAILABLE"         "seq_hash512 +UNAVAILABLE" bash "$QC/verify_reference_fingerprint.sh" "$W/info_na.json" "$W/fp3.json"
bad "wrong requested k detected"                "MISMATCH \(requested 25\)" "$PY" "$TOOLS" verify-fingerprint --observed "$INDEX_DIR" --expect-k 25 --require seq_hash --expect-seq-hash "$EXPECT_SEQ_HASH"

echo "== specimen eligibility (validated_pairs.tsv + exclusions)"
bad "11-37 blocked pending provenance adjudication" "excluded: duplicate Drive objects" env SPECIMEN=GM37 bash "$Q"
bad "11-37 blocked by short id too"             "excluded: duplicate" env SPECIMEN=11-37 bash "$Q"
bad "lift without reference rejected"           "needs CORE_ID=<adjudication reference>" env SPECIMEN=GM37 LIFT_EXCLUSIONS="20088XR-11-37=" bash "$Q"
bad "corrupt 11-43 blocked"                     "excluded: corrupt FASTQ" env SPECIMEN=GM43 bash "$Q"
bad "specimen absent from validated_pairs blocked" "only gzip-validated pairs are eligible" env SPECIMEN=GM03 bash "$Q"
bad "extra EXCLUDED_SPECIMENS honoured"         "excluded by EXCLUDED_SPECIMENS" env SPECIMEN=GM01 EXCLUDED_SPECIMENS="GM01" bash "$Q"
mkdir -p "$W/vp"; sed 's/^20088XR-11-01/20088XR-11-09/' "$F/validated_pairs.tsv" > "$W/vp/label.tsv"
bad "inconsistent validated label blocked"      "inconsistent with" env SPECIMEN=GM01 VALIDATED_PAIRS="$W/vp/label.tsv" bash "$Q"
mkdir -p "$W/rawcopy"; cp "$F/raw/20088XR-11-01_S0_L001_R1_001.fastq.gz" "$F/raw/20088XR-11-01_S0_L001_R2_001.fastq.gz" "$W/rawcopy/"
# both mates move together: read_validated_pairs requires R1 and R2 in the same directory
sed "s#$F/raw/20088XR-11-01_S0_L001_R\([12]\)_001#$W/rawcopy/20088XR-11-01_S0_L001_R\1_001#g" "$F/validated_pairs.tsv" > "$W/vp/otherfile.tsv"
bad "validated path pointing at a different file blocked" "are different files" env SPECIMEN=GM01 VALIDATED_PAIRS="$W/vp/otherfile.tsv" bash "$Q"
printf 'sample\tr1\tr2\n20088XR-11-01\tbroken\n' > "$W/vp/malformed.tsv"
bad "malformed validated_pairs rejected"        "malformed line|no pairs" env SPECIMEN=GM01 VALIDATED_PAIRS="$W/vp/malformed.tsv" bash "$Q"
sed 's/^GM02,/..,/' "$F/samples.csv" > "$W/vp/dotdot.csv"
bad "traversal-like sample_id rejected"         "unsafe sample_id" "$PY" "$TOOLS" resolve-specimen --sheet "$W/vp/dotdot.csv" --validated-pairs "$F/validated_pairs.tsv" --specimen 20088XR-11-02
sed 's#20088XR-11-02_S0_L001_R1_001.fastq.gz,#../x/20088XR-11-02_S0_L001_R1_001.fastq.gz,#' "$F/samples.csv" > "$W/vp/travfile.csv"
bad "traversal in raw filename rejected"        "unsafe r1" "$PY" "$TOOLS" resolve-specimen --sheet "$W/vp/travfile.csv" --validated-pairs "$F/validated_pairs.tsv" --specimen GM02

echo "== r3 manifest parser: 48-specimen eligibility replica (real core-ID set and historical layout)"
# Same specimen set as the committed sample sheet: 11-01..11-44, 11-46..11-48, 11-50 (48); validated pairs exclude the
# corrupt 11-43/44/46/47/48 (43 rows). Expected: 42 eligible; GM37 + the five corrupt specimens blocked.
E="$W/elig48"; mkdir -p "$E"; ERAW=/data/project/raw_fastq
ECORES=""; for n in $(seq 1 44) 46 47 48 50; do ECORES="$ECORES $(printf '%02d' "$n")"; done
{ echo "sample_id,core_id,project,species,tissue,specimen,disease_label,lane,index_s,r1,r2,md5_r1,md5_r2"
  for n in $ECORES; do c=20088XR-11-$n; echo "GM$n,$c,20088-11,human,melanoma,FFPE,melanoma,L004,S1$n,${c}_S1${n}_L004_R1_001.fastq.gz,${c}_S1${n}_L004_R2_001.fastq.gz,$(printf '%032d' "1$n"),$(printf '%032d' "2$n")"; done; } > "$E/sheet.csv"
evp() {  # evp FILE LABELSTYLE(core|short|gm) [noheader]: write a 43-row validated_pairs replica
  { [ "${3:-}" = noheader ] || printf 'sample\tr1_path\tr2_path\n'
    for n in $ECORES; do case $n in 43|44|46|47|48) continue ;; esac; c=20088XR-11-$n
      case $2 in core) l=$c ;; short) l=11-$n ;; gm) l=GM$n ;; esac
      printf '%s\t%s\t%s\n' "$l" "$ERAW/${c}_S1${n}_L004_R1_001.fastq.gz" "$ERAW/${c}_S1${n}_L004_R2_001.fastq.gz"; done; } > "$1"; }
eligible() {  # eligible VP: print "<n_eligible> <blocked ids>" over all 48 sheet specimens
  local n=0 b="" s
  for s in $ECORES; do if "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$1" --specimen "GM$s" >/dev/null 2>"$E/err.txt"
    then n=$((n+1)); else b="$b GM$s"; fi; done; echo "$n$b"; }
evp "$E/core.tsv" core; evp "$E/short.tsv" short; evp "$E/gm.tsv" gm; evp "$E/nohdr.tsv" core noheader
for v in core short gm nohdr; do
  cond "historical layout ($v labels/header): 42 eligible, GM37 + corrupt blocked" test "$(eligible "$E/$v.tsv")" = "42 GM37 GM43 GM44 GM46 GM47 GM48"; done
bad "GM37 blocked by exclusion, not by manifest"   "excluded: duplicate Drive objects" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/core.tsv" --specimen GM37
bad "corrupt GM44 blocked"                          "excluded: corrupt FASTQ" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/core.tsv" --specimen GM44
bad "lifted corrupt GM44 still not eligible (not a validated pair)" "only gzip-validated pairs are eligible" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/core.tsv" --specimen GM44 --lift-exclusion "20088XR-11-44=TEST"
"$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/core.tsv" --specimen 11-02 > "$E/res02.txt"
has "resolved output unchanged in shape (8 fields; basenames, sheet MD5s, label)" $'^GM02\t20088XR-11-02\t20088XR-11-02_S102_L004_R1_001.fastq.gz\t20088XR-11-02_S102_L004_R2_001.fastq.gz\t0+102\t0+202\t20088XR-11-02\tnone$' "$E/res02.txt"
# manifest-wide fail-closed: any defective row blocks every specimen (documented r3 behaviour change from r2)
sed '21s/^/ /' "$E/core.tsv" > "$E/padded.tsv"
cond "one padded label: whole manifest refused (0 eligible)" test "$(eligible "$E/padded.tsv" | cut -d' ' -f1)" = 0
bad "  ... reason reported"                          "empty field or surrounding whitespace" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/padded.tsv" --specimen GM01
{ cat "$E/core.tsv"; grep '20088XR-11-37_' "$E/core.tsv"; } > "$E/dup37.tsv"
bad "11-37 listed twice: duplicate refused manifest-wide" "duplicate (specimen|sample label)" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/dup37.tsv" --specimen GM01
sed "s#$ERAW/##g" "$E/core.tsv" > "$E/relative.tsv"
bad "relative manifest paths refused"               "path is not absolute" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/relative.tsv" --specimen GM01
sed '1s/^sample/id/' "$E/core.tsv" > "$E/badheader.tsv"
bad "unknown header word refused"                   "malformed line 1" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/badheader.tsv" --specimen GM01
awk -F'\t' 'NR==3 {t=$2; $2=$3; $3=t} {print}' OFS='\t' "$E/core.tsv" > "$E/swapped.tsv"
bad "R1/R2 swapped refused"                         "R2 is not the mate of R1|does not end in _R1_" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/swapped.tsv" --specimen GM01
printf 'sample\tr1_path\tr2_path\r\n' > "$E/crlf.tsv"; tail -n +2 "$E/core.tsv" | sed 's/$/\r/' >> "$E/crlf.tsv"
bad "CRLF manifest refused"                         "CR \(Windows\)" "$PY" "$TOOLS" resolve-specimen --sheet "$E/sheet.csv" --validated-pairs "$E/crlf.tsv" --specimen GM01

echo "== quantify_specimen.sbatch: path, lock and adoption safety"
sed 's/,\([0-9a-f]\{32\}\)$/,ffffffffffffffffffffffffffffffff/' "$F/samples.csv" > "$W/vp/bad_md5.csv"
bad "quant fails on raw MD5 mismatch"           "MD5 mismatch for .*R2_001" env SPECIMEN=GM01 SAMPLE_SHEET="$W/vp/bad_md5.csv" bash "$Q"
cond "MD5 failure leaves no specimen output"    test ! -e "$QUANT_ROOT/GM01"
cond "lock released after failure"              test ! -e "$QUANT_ROOT/.GM01.lock"
bad "quant fails on container SHA mismatch"     "SHA-256 mismatch for .*fastp.sif" env SPECIMEN=GM01 FASTP_SIF_SHA256=$(printf '%064d' 0) bash "$Q"
mkdir -p "$W/raw_rw/sub"; cp "$F"/raw/* "$W/raw_rw/"
sed "s#$F/raw#$W/raw_rw#g" "$F/validated_pairs.tsv" > "$W/vp/rw.tsv"
bad "quant refuses QUANT_ROOT inside RAW_DIR"   "must not be inside RAW_DIR" env SPECIMEN=GM01 RAW_DIR="$W/raw_rw" VALIDATED_PAIRS="$W/vp/rw.tsv" QUANT_ROOT="$W/raw_rw/sub" bash "$Q"
ln -s "$W/raw_rw/sub" "$W/qlink"
bad "quant refuses symlinked QUANT_ROOT"        "not a real directory" env SPECIMEN=GM01 RAW_DIR="$W/raw_rw" VALIDATED_PAIRS="$W/vp/rw.tsv" QUANT_ROOT="$W/qlink" bash "$Q"
bad "quant refuses RAW_DIR inside QUANT_ROOT"   "RAW_DIR must not be inside QUANT_ROOT" env SPECIMEN=GM01 RAW_DIR="$W/raw_rw" VALIDATED_PAIRS="$W/vp/rw.tsv" QUANT_ROOT="$W" bash "$Q"
ln -s "$F/raw" "$W/rawlink"
bad "quant refuses symlinked RAW_DIR"           "not a real directory" env SPECIMEN=GM01 RAW_DIR="$W/rawlink" bash "$Q"
mkdir -p "$W/ref/idx.partial.9"; cp "$INDEX_DIR"/*.json "$W/ref/idx.partial.9/"
bad "quant refuses a partial index"             "partial \(unadopted\) index" env SPECIMEN=GM01 INDEX_DIR="$W/ref/idx.partial.9" bash "$Q"
mkdir -p "$W/ref/idx_bad_fp"; cp "$W/info_mm.json" "$W/ref/idx_bad_fp/info.json"
bad "quant refuses index with wrong fingerprint" "index directory fingerprints did not all MATCH" env SPECIMEN=GM01 INDEX_DIR="$W/ref/idx_bad_fp" bash "$Q"
rm -rf "$QUANT_ROOT"/.GM01.partial.*
mkdir -p "$QUANT_ROOT/.GM01.lock"
bad "concurrent job blocked by specimen lock"   "lock held or stale" env SPECIMEN=GM01 bash "$Q"
cond "lock held by other job is not removed"    test -d "$QUANT_ROOT/.GM01.lock"
rmdir "$QUANT_ROOT/.GM01.lock"
bad "truncated quant.sf never adopted"          "failed structural validation" env SPECIMEN=GM01 STUB_TRUNCATE=1 bash "$Q"
has "  ... reason recorded"                     "PROBLEM quant.sf has 99 rows" "$(ls "$QUANT_ROOT"/.GM01.partial.*/structure_check.txt | head -1)"
rm -rf "$QUANT_ROOT"/.GM01.partial.*
bad "wrong library orientation never adopted"   "failed structural validation" env SPECIMEN=GM01 STUB_LIBTYPE=ISF bash "$Q"
cond "failed runs left no adopted output"       test ! -e "$QUANT_ROOT/GM01"
rm -rf "$QUANT_ROOT"/.GM01.partial.*

echo "== quantify_specimen.sbatch: actionable diagnostics (no silent exits)"
bad "unloadable APPTAINER_MODULES reported with Lmod message" "module load failed for step modules 'nonexistent/9'.*Lmod has detected" \
    env SPECIMEN=GM01 RQC_SKIP_MODULE_LOAD=0 APPTAINER_MODULES=nonexistent/9 bash "$Q"
bad "unloadable SALMON_MODULES reported"        "module load failed for step modules 'nonexistent/8'" \
    env SPECIMEN=GM01 RQC_SKIP_MODULE_LOAD=0 SALMON_MODULES=nonexistent/8 bash "$Q"
bad "foreign Lmod (BASH_ENV) failure is diagnosed, not silent" "fastp version check failed \(exit 255\).*image format not recognized" \
    env SPECIMEN=GM01 RQC_SKIP_MODULE_LOAD=0 BASH_ENV="$FOREIGN/init.sh" bash "$Q"
has "  ... diagnostic names the resolved executable" "resolved=$FOREIGN/bin/apptainer" "$T/last.log"
has "  ... the foreign module really was the one called" '^FOREIGN load apptainer$' "$FOREIGN/calls.log"
cond "diagnosed failures left no output or lock" test ! -e "$QUANT_ROOT/GM01" -a ! -e "$QUANT_ROOT/.GM01.lock"
rm -rf "$QUANT_ROOT"/.GM01.partial.*

echo "== quantify_specimen.sbatch: success, isolation and provenance"
ok  "quant GM01 (Lmod stub, isolated modules)" env SPECIMEN=GM01 RQC_SKIP_MODULE_LOAD=0 bash "$Q"
has "  ... stub module loaded apptainer for fastp" '^load apptainer$' "$RQC_STUB_MODULE_LOG"
P1="$QUANT_ROOT/GM01/provenance.json"; P1F=$( { [ -f "$P1" ] && cp "$P1" "$W/p1.json"; } ; flat "$W/p1.json")
has "salmon step ran with bzip2/1.0.6"          '"salmon_step_modules": "[^"]*bzip2/1\.0\.6' "$P1"
hasnt "salmon step never saw bzip2/1.0.8"       '"salmon_step_modules": "[^"]*bzip2/1\.0\.8' "$P1"
has "fastp step recorded with apptainer modules" '"fastp_step_modules": "[^"]*apptainer/1\.4\.5' "$P1"
ok  "quant GM02 by short validated label"      env SPECIMEN=11-02 bash "$Q"
ok  "quant 11-37 only with recorded lift"      env SPECIMEN=GM37 LIFT_EXCLUSIONS="20088XR-11-37=ADJ-2026-001" bash "$Q"
has "  ... lift recorded in provenance"         '"lifted_exclusions": "20088XR-11-37=ADJ-2026-001"' "$QUANT_ROOT/GM37/provenance.json"
cond "outputs present, trimmed reads deleted"   test -s "$QUANT_ROOT/GM01/quant/quant.sf" -a -z "$(find "$QUANT_ROOT/GM01" -name '*.fq.gz')"
cond "no lock or partial left after success"    test -z "$(ls -A "$QUANT_ROOT" | grep -E '^\.')"
has "provenance status complete"                '"status": "complete"' "$P1"
has "fastp version labelled candidate"          '"fastp_version": \{[^}]*"status": "candidate"' "$P1F"
has "salmon args labelled verified-collaborator" '"salmon_quant_args": \{[^}]*"status": "verified-collaborator"' "$P1F"
has "fastp threads labelled inferred"           '"fastp_threads": \{[^}]*"status": "inferred"' "$P1F"
has "index k observed-this-run"                 '"index_kmer": \{[^}]*"status": "observed-this-run"' "$P1F"
has "  ... value 31, collaborator k unverified" '"index_kmer": \{[^}]*remains unverified[^}]*"value": "31"' "$P1F"
bad "quant refuses output collision"            "refusing to overwrite.*GM01" env SPECIMEN=GM01 bash "$Q"
ok  "validate-quant-dir accepts adopted output" "$PY" "$TOOLS" validate-quant-dir --dir "$QUANT_ROOT/GM01" --specimen GM01 --expect-targets 2000 --expect-libtype ISR
RAW_AFTER=$(for f in "$F"/raw/*; do h256 "$f"; done | tr '\n' ' ')
cond "raw inputs unchanged"                     test "$RAW_BEFORE" = "$RAW_AFTER"

echo "== compare-quant (rule provisional-v1)"
R="$W/refset"; mkdir -p "$R"
cp "$F/reference_quant.sf" "$R/quant.sf"; gzip -c "$R/quant.sf" > "$R/quant.sf.gz"
cp "$QUANT_ROOT/GM01/quant/aux_info/meta_info.json" "$R/meta_info.json"; cp "$QUANT_ROOT/GM01/quant/lib_format_counts.json" "$R/lib_format_counts.json"
chmod u+w "$R"/*
RAWN=$("$PY" -c "import json,sys;print(json.load(open(sys.argv[1]))['summary']['before_filtering']['total_reads'])" "$QUANT_ROOT/GM01/fastp/GM01.json")
printf 'sample_id\tcore_id\treads_raw\treads_after_trim\nGM01\t20088XR-11-01\t%s\t%s\nGM02\t20088XR-11-02\t%s\t%s\n' "$RAWN" "$RAWN" "$RAWN" "$RAWN" > "$R/run_qc.tsv"
cmp() {  # cmp [VAR=value ...] [extra args]: VARs override CS CT CQ CM CL CR
  local CS=GM01 CT="$QUANT_ROOT/GM01" CQ="$R/quant.sf.gz" CM="$R/meta_info.json" CL="$R/lib_format_counts.json" CR="$R/run_qc.tsv"
  while [ $# -gt 0 ] && [[ "$1" =~ ^C[STQMLR]= ]]; do local "$1"; shift; done
  "$PY" "$TOOLS" compare-quant --specimen "$CS" --test-dir "$CT" --reference-quant "$CQ" --reference-meta "$CM" \
    --reference-lib-format "$CL" --reference-run-qc "$CR" --gene-map-fasta "$F/genemap.fa.gz" --expect-targets 2000 "$@"; }
ok  "identical reproduction: PASS"             cmp --json-out "$W/cmp1.json"
has "report: rule provisional-v1"               '"version": "provisional-v1"' "$W/cmp1.json"
has "report: gene level separate, not gated"    '"gene_level_reported_not_gated"' "$W/cmp1.json"
has "report: transcript level separate"         '"transcript_level"' "$W/cmp1.json"
has "report: mapped fragments per raw pair"     'mapped_fragments_per_raw_pair_test' "$W/cmp1.json"
bad "wrong specimen identity: FAIL"             "gate specimen_identity +FAIL" cmp CS=GM02
"$PY" -c "import sys; r=open(sys.argv[1]).read().replace('\t$RAWN\t$RAWN\n', '\t999\t$RAWN\n', 1); open(sys.argv[2],'w').write(r)" "$R/run_qc.tsv" "$W/run_qc_raw.tsv"
bad "raw input read count differs: FAIL"        "gate raw_input_reads_identical +FAIL" cmp CR="$W/run_qc_raw.tsv"
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); d['expected_format']='ISF'; json.dump(d,open(sys.argv[2],'w'))" "$R/lib_format_counts.json" "$W/lib_isf.json"
bad "library orientation differs: FAIL"         "gate library_format_matches +FAIL" cmp CL="$W/lib_isf.json"
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); d['num_valid_targets']=1999; json.dump(d,open(sys.argv[2],'w'))" "$R/meta_info.json" "$W/meta_targets.json"
bad "target count differs: FAIL"                "gate target_counts_equal +FAIL" cmp CM="$W/meta_targets.json"
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); d['index_name_hash512']='00'; json.dump(d,open(sys.argv[2],'w'))" "$R/meta_info.json" "$W/meta_fp.json"
bad "reference SHA-512 fingerprint differs: FAIL" "gate index_name_hash512_identical +FAIL" cmp CM="$W/meta_fp.json"
"$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); d['percent_mapped']+=2.5; json.dump(d,open(sys.argv[2],'w'))" "$R/meta_info.json" "$W/meta_shift.json"
bad "mapping-rate shift > 1 pp: FAIL"           "gate abs_mapping_diff_le_1pp +FAIL" cmp CM="$W/meta_shift.json"
sed '1500,$d' "$F/reference_quant.sf" > "$W/truncref.sf"
bad "transcript ID set differs: FAIL"           "gate identical_transcript_ids +FAIL" cmp CQ="$W/truncref.sf"
mkdir -p "$W/alt"; cp -R "$QUANT_ROOT/GM01" "$W/alt/GM01"; chmod -R u+w "$W/alt"
"$PY" - "$W/alt/GM01/quant/quant.sf" <<'P'
import random,sys
random.seed(7); L=open(sys.argv[1]).read().splitlines()
out=[L[0]]+["\t".join(f[:3]+[f"{float(f[3])*random.choice([0.05,1,20]):.6f}",f[4]]) for f in (l.split("\t") for l in L[1:])]
open(sys.argv[1],"w").write("\n".join(out)+"\n")
P
bad "perturbed expression: FAIL"                "gate transcript_pearson_log1p_tpm_ge_0.99 +FAIL" cmp CT="$W/alt/GM01"
cp -R "$QUANT_ROOT/GM01" "$W/alt/.GM01.partial.77"
bad "partial work dir as test: FAIL"            "PROBLEM path is a partial work directory" cmp CT="$W/alt/.GM01.partial.77"
cp -R "$QUANT_ROOT/GM01" "$W/alt/GM01_incomplete"; chmod -R u+w "$W/alt"
sed 's/"status": "complete"/"status": "running"/' "$W/alt/GM01_incomplete/provenance.json" > "$W/alt/p.json" && mv "$W/alt/p.json" "$W/alt/GM01_incomplete/provenance.json"
bad "incomplete provenance: FAIL"               "PROBLEM provenance status is 'running'" cmp CT="$W/alt/GM01_incomplete"
bad "threshold change without new rule version rejected" "require --rule-version" cmp --min-pearson 0.95
bad "threshold change without justification rejected"    "require --rule-version" cmp --min-pearson 0.95 --rule-version alt-v1
ok  "versioned, justified alternative rule runs" cmp --min-pearson 0.95 --rule-version alt-v1 --rule-justification "test" --json-out "$W/cmp_alt.json"
has "  ... alternative rule labelled in report" '"version": "alt-v1"' "$W/cmp_alt.json"
bad "compare refuses report collision"          "refusing to overwrite" cmp --json-out "$W/cmp1.json"

echo "== audit-historical (pre-provenance-schema run; read-only)"
HD="$W/hist/GM01_candidate_99"; mkdir -p "$HD/aux_info" "$HD/logs"
cp "$QUANT_ROOT/GM01/quant/quant.sf" "$HD/quant.sf"; cp "$QUANT_ROOT/GM01/quant/aux_info/meta_info.json" "$HD/aux_info/"
cp "$QUANT_ROOT/GM01/quant/lib_format_counts.json" "$HD/"; cp "$QUANT_ROOT/GM01/fastp/GM01.json" "$HD/logs/GM01.fastp.json"
printf '{"salmon_version":"1.10.3","libType":"A","mates1":"/x/tmp/GM01_R1.fq.gz","mates2":"/x/tmp/GM01_R2.fq.gz","threads":"6","validateMappings":[],"gcBias":[],"seqBias":[]}\n' > "$HD/cmd_info.json"
chmod -R a-w "$HD"
HD_BEFORE=$(cd "$HD" && find . -type f | sort | while read -r f; do h256 "$f"; done | tr '\n' ' ')
(cd "$R" && for f in quant.sf.gz meta_info.json lib_format_counts.json run_qc.tsv; do printf '%s  %s\n' "$(h256 "$f")" "$f"; done) > "$R/MANIFEST.sha256"
hist() {  # hist [HD=dir] [extra args]
  local D="$HD"; [[ "${1:-}" == HD=* ]] && { D=${1#HD=}; shift; }
  "$PY" "$TOOLS" audit-historical --specimen GM01 --test-dir "$D" --reference-quant "$R/quant.sf.gz" --reference-meta "$R/meta_info.json" \
    --reference-lib-format "$R/lib_format_counts.json" --reference-run-qc "$R/run_qc.tsv" --gene-map-fasta "$F/genemap.fa.gz" \
    --expect-targets 2000 "$@"; }
bad "validate-quant-dir still rejects the ad hoc historical layout" "^INVALID$" "$PY" "$TOOLS" validate-quant-dir --dir "$HD" --specimen GM01
mkdir -p "$W/hist/stdlayout"; cp -R "$QUANT_ROOT/GM01/quant" "$QUANT_ROOT/GM01/fastp" "$W/hist/stdlayout/"; chmod -R u+w "$W/hist/stdlayout"
bad "validate-quant-dir still rejects a run without provenance" "PROBLEM missing provenance.json" "$PY" "$TOOLS" validate-quant-dir --dir "$W/hist/stdlayout" --specimen GM01
ok  "audit discovers the standard layout too (no cmd_info: identity UNAVAILABLE -> exit 2 tolerated below)" \
    bash -c '"$@"; [ $? -eq 2 ]' _ "$PY" "$TOOLS" audit-historical --specimen GM01 --test-dir "$W/hist/stdlayout" --reference-quant "$R/quant.sf.gz" \
    --reference-meta "$R/meta_info.json" --reference-lib-format "$R/lib_format_counts.json" --reference-run-qc "$R/run_qc.tsv" \
    --gene-map-fasta "$F/genemap.fa.gz" --expect-targets 2000
has "  ... all provisional-v1 gates PASS, identity UNAVAILABLE" "identity +specimen_identity_from_recorded_inputs +UNAVAILABLE" "$T/last.log"
hasnt "  ... and no gate FAIL"                     "^gate .*FAIL$" "$T/last.log"
ok  "historical concordant run: HISTORICAL_CONCORDANT" hist --json-out "$W/hist1.json"
has "  ... provenance reported absent, not fabricated" '"status": "absent: run predates the provenance schema; not fabricated"' "$W/hist1.json"
has "  ... no adoption"                           '"adoption": "none: historical audit only' "$W/hist1.json"
has "  ... rule provisional-v1"                   '"version": "provisional-v1"' "$W/hist1.json"
has "  ... identity from recorded inputs PASS"    '"specimen_identity_from_recorded_inputs": "PASS"' "$W/hist1.json"
has "  ... salmon/fastp consistency PASS"         '"salmon_processed_equals_fastp_retained_pairs": "PASS"' "$W/hist1.json"
has "  ... reference manifest verified"           '"reference_bundle_manifest_verified": "PASS"' "$W/hist1.json"
has "  ... gene level separate, not gated"        '"gene_level_reported_not_gated"' "$W/hist1.json"
has "  ... mapping yield per raw pair"            'mapped_fragments_per_raw_pair_test' "$W/hist1.json"
has "  ... evidence files hashed"                 '"evidence_files"' "$W/hist1.json"
HD_AFTER=$(cd "$HD" && find . -type f | sort | while read -r f; do h256 "$f"; done | tr '\n' ' ')
cond "audited directory unchanged (read-only)"    test "$HD_BEFORE" = "$HD_AFTER" -a "$(find "$HD" -type f | wc -l | tr -d ' ')" = 5
bad "json-out inside audited dir refused"         "must be outside the audited directory" hist --json-out "$HD/logs/x.json"
bad "historical report collision refused"         "refusing to overwrite" hist --json-out "$W/hist1.json"
bad "threshold change without new rule version rejected" "require --rule-version" hist --min-pearson 0.9
mkdir -p "$W/hist/nofastp/aux_info"; cp "$HD/quant.sf" "$HD/lib_format_counts.json" "$HD/cmd_info.json" "$W/hist/nofastp/"; cp "$HD/aux_info/meta_info.json" "$W/hist/nofastp/aux_info/"
bad "missing fastp JSON: HISTORICAL_INCONCLUSIVE (exit 2)" "raw_input_reads_identical +UNAVAILABLE" hist HD="$W/hist/nofastp"
has "  ... verdict inconclusive, not concordant"   "VERDICT HISTORICAL_INCONCLUSIVE" "$T/last.log"
cp -R "$W/hist/nofastp" "$W/hist/wrongspec"; chmod -R u+w "$W/hist/wrongspec"; cp "$HD/logs/GM01.fastp.json" "$W/hist/wrongspec/"
sed 's/GM01_R/GM02_R/g' "$HD/cmd_info.json" > "$W/hist/wrongspec/cmd_info.json"
bad "inputs recorded for another specimen: DISCORDANT" "specimen_identity_from_recorded_inputs +FAIL" hist HD="$W/hist/wrongspec"
cp -R "$W/hist/wrongspec" "$W/hist/twoquant"; chmod -R u+w "$W/hist/twoquant"; cp "$HD/cmd_info.json" "$W/hist/twoquant/"
mkdir "$W/hist/twoquant/old"; cp "$HD/quant.sf" "$W/hist/twoquant/old/"
bad "two candidate quant.sf files: ambiguous, refused" "ambiguous: 2 candidate quant files" hist HD="$W/hist/twoquant"
ok  "  ... explicit --test-quant resolves it"      hist HD="$W/hist/twoquant" --test-quant "$W/hist/twoquant/quant.sf"
mkdir -p "$W/hist/nometa"; cp "$HD/quant.sf" "$HD/lib_format_counts.json" "$W/hist/nometa/"
bad "missing meta_info.json refused"              "required meta file not found" hist HD="$W/hist/nometa"
cp -R "$W/hist/twoquant" "$W/hist/perturbed"; chmod -R u+w "$W/hist/perturbed"; rm -rf "$W/hist/perturbed/old"; cp "$W/alt/GM01/quant/quant.sf" "$W/hist/perturbed/quant.sf"
bad "perturbed expression: DISCORDANT"            "gate +transcript_pearson_log1p_tpm_ge_0.99 +FAIL" hist HD="$W/hist/perturbed"
cp -R "$R" "$W/refbad"; chmod -R u+w "$W/refbad"; printf 'x\n' >> "$W/refbad/run_qc.tsv"
bad "reference bundle not matching its manifest: FAIL" "reference_bundle_manifest_verified +FAIL" hist --reference-run-qc "$W/refbad/run_qc.tsv" --reference-manifest "$R/MANIFEST.sha256"
bad "wrong gene-map MD5 refused"                  "gene map FASTA MD5" hist --gene-map-md5 00000000000000000000000000000000

echo
echo "local tests: $pass passed, $fail failed  (stubs only; not evidence of Rivanna behaviour)"
[ "$fail" -eq 0 ]
