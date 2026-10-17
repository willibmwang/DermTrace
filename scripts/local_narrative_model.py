"""Run the frozen narrative prompt on a pinned local Transformers model."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path


def generate_local_attempts(records, repeats, model_path, batch_size=16):
    import torch
    import transformers
    from transformers import AutoModelForCausalLM, AutoTokenizer
    import demo_utils as du

    if not torch.cuda.is_available():
        raise RuntimeError('Local narrative audit requires an allocated CUDA GPU.')
    model_path = Path(model_path).resolve()
    if not model_path.is_dir():
        raise ValueError('Use an existing, pinned model snapshot directory.')
    tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True, padding_side='left')
    if tokenizer.pad_token_id is None:
        tokenizer.pad_token_id = tokenizer.eos_token_id
    model = AutoModelForCausalLM.from_pretrained(
        model_path, local_files_only=True, torch_dtype=torch.bfloat16,
        attn_implementation='sdpa',
    ).to('cuda').eval()
    model_id = 'Qwen/Qwen3-1.7B@' + model_path.name
    metadata = {
        'model': model_id, 'snapshot': model_path.name,
        'config_sha256': hashlib.sha256((model_path / 'config.json').read_bytes()).hexdigest(),
        'transformers': transformers.__version__, 'torch': torch.__version__,
        'device': torch.cuda.get_device_name(), 'dtype': 'bfloat16',
        'enable_thinking': False, 'do_sample': False, 'max_new_tokens': 400,
        'batch_size': batch_size,
    }
    tasks = [(record, repeat) for record in records for repeat in range(repeats)]
    for offset in range(0, len(tasks), batch_size):
        batch = tasks[offset:offset + batch_size]
        payloads = []
        for record, repeat in batch:
            packet = du.external_evidence_packet(record['packet'])
            payloads.append({
                'model': model_id,
                'messages': [{'role': 'user', 'content': du._copilot_prompt(packet)}],
                'temperature': 0, 'max_tokens': 400,
            })
        rendered = [tokenizer.apply_chat_template(
            payload['messages'], tokenize=False, add_generation_prompt=True,
            enable_thinking=False,
        ) for payload in payloads]
        encoded = tokenizer(rendered, return_tensors='pt', padding=True).to('cuda')
        started = time.monotonic()
        try:
            with torch.inference_mode():
                generated = model.generate(
                    **encoded, max_new_tokens=400, do_sample=False,
                    pad_token_id=tokenizer.pad_token_id,
                )
            new_tokens = generated[:, encoded.input_ids.shape[1]:].cpu()
            texts = tokenizer.batch_decode(new_tokens, skip_special_tokens=True)
            for index, ((record, repeat), payload, narrative) in enumerate(zip(batch, payloads, texts)):
                token_ids = new_tokens[index].tolist()
                ended = tokenizer.eos_token_id in token_ids
                response = {
                    **metadata, 'id': f'local-{offset + index:04d}',
                    'finish_reason': 'stop' if ended else 'length',
                    'rendered_prompt_sha256': hashlib.sha256(rendered[index].encode()).hexdigest(),
                    'usage': {'prompt_tokens': int(encoded.attention_mask[index].sum()),
                              'completion_tokens': token_ids.index(tokenizer.eos_token_id) + 1 if ended else len(token_ids)},
                }
                yield narrative.strip(), [payload, {'response_metadata': response}], None
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            for payload in payloads:
                yield '', [payload, {'response_metadata': metadata}], 'OutOfMemoryError'
        print(f'Generated batch {offset + 1}-{min(offset + batch_size, len(tasks))}/{len(tasks)} '
              f'in {time.monotonic() - started:.1f}s', flush=True)

