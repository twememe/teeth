#!/usr/bin/env python3
"""Prepare isolated Full181 inputs and one matched A/H/L ColabFold batch.

No old inputs are written, and no hashes or model smoke tests are performed.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/full181_supplement'
LOG = ROOT / 'logs/full181_supplement'
EXPECTED = ('MSWGLQILPCLSLILLLWNQVPGLEGQEFRSGSCQVTGVVLPELWEAFWTVKNTVQTQDD'
            'ITSIRLLKPQVLRNVSGAESCYLAHSLLKFYLNTVFKNYHSKIAKFKVLRSFSTLANNFI'
            'VIMSQLQPSKDNSMLPISESAHQRFLLFRRAFKQLDTEVALVKAFGEVDILLTWMQKFYHL')


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                     prefix=path.name + '.tmp_', delete=False) as handle:
        handle.write(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
        temporary = Path(handle.name)
    temporary.replace(path)


def metadata(path):
    path = Path(path).resolve()
    if not path.is_file():
        return {'path': str(path), 'exists': False}
    stat = path.stat()
    return {'path': str(path), 'exists': True, 'bytes': stat.st_size,
            'mtime_ns': stat.st_mtime_ns, 'fingerprint_method': 'path_size_mtime_no_hash'}


def fasta(path):
    records = {}
    key = None
    for raw in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line = raw.strip()
        if line.startswith('>'):
            key = line[1:].split()[0]
            if key in records:
                raise ValueError(f'duplicate FASTA record: {key}')
            records[key] = ''
        elif line:
            if key is None:
                raise ValueError('FASTA sequence before header')
            records[key] += line.upper()
    return records


def input_sequences():
    full = next(iter(fasta(ROOT / 'data/IL24_full_1_181.fasta').values()))
    if full != EXPECTED:
        differences = [{'position': i + 1, 'expected': a, 'actual': b}
                       for i, (a, b) in enumerate(__import__('itertools').zip_longest(EXPECTED, full)) if a != b]
        write_json(OUT / 'input_sequence_mismatch.json', differences)
        raise ValueError('Full sequence mismatch; affected GPU tasks stopped')
    ab = fasta(ROOT / 'data/VH_VL_clean.fasta')
    vh = next(v for k, v in ab.items() if 'VH' in k)
    vl = next(v for k, v in ab.items() if 'VL' in k)
    if (len(full), len(vh), len(vl)) != (181, 119, 106):
        raise ValueError('A/H/L length contract violated')
    return dict(zip('AHL', (full, vh, vl)))


def prepare():
    seqs = input_sequences()
    inp = OUT / 'inputs'
    inp.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    if (inp / 'sequences.json').exists() and json.loads((inp / 'sequences.json').read_text(encoding='utf-8')) != seqs:
        raise ValueError('Existing supplement input sequences differ; refusing overwrite')
    write_json(inp / 'sequences.json', seqs)
    for condition, target in [('empty', 'full_blind'), ('formal_msa', 'full_blind_msa')]:
        text = 'version: 1\nsequences:\n'
        for chain, sequence in seqs.items():
            msa = 'empty' if condition == 'empty' else f'results/full181_supplement/msa/{chain}.csv'
            text += f'  - protein:\n      id: {chain}\n      sequence: {sequence}\n      msa: {msa}\n'
        yaml_path = inp / f'{target}.yaml'
        if yaml_path.exists() and yaml_path.read_text(encoding='utf-8') != text:
            raise ValueError(f'Existing input differs; refusing overwrite: {yaml_path}')
        if not yaml_path.exists():
            yaml_path.write_text(text, encoding='utf-8')
    params = {'model': 'boltz1', 'boltz_version': '2.2.1', 'accelerator': 'gpu',
              'devices': 1, 'diffusion_samples': 1, 'recycling_steps': 3,
              'sampling_steps': 200, 'no_kernels': True, 'num_workers': 0,
              'output_format': 'pdb', 'write_full_pae': True, 'seeds': list(range(1, 16)),
              'hash_policy': 'not_computed_per_direct_user_instruction',
              'resume_policy': 'exact input text and file metadata plus structure/confidence/PAE validation'}
    write_json(OUT / 'parameters.json', params)
    if not (OUT / 'baseline_file_metadata.json').exists():
        old = [metadata(p) for folder in ('data', 'results', 'src', 'docs', 'models', 'logs')
               for p in (ROOT / folder).rglob('*') if p.is_file()
               and 'full181_supplement' not in p.parts and '__pycache__' not in p.parts
               and not p.name.startswith(('full181_', 'prepare_full181_', 'run_full181_', 'analyze_full181_', 'summarize_full181_'))]
        write_json(OUT / 'baseline_file_metadata.json', old)
    tasks = []
    for condition, target in [('empty', 'full_blind'), ('formal_msa', 'full_blind_msa')]:
        for seed in range(1, 16):
            tasks.append({'construct': 'Full181', 'condition': condition, 'seed': seed,
                          'task_id': f'Full181_{condition}_seed{seed}',
                          'source_path': f'results/full181_supplement/raw/{condition}/seed{seed}',
                          'input_path': f'results/full181_supplement/inputs/{target}.yaml',
                          'status': 'planned', 'failure_reason': '',
                          'input_hash': 'not_computed_user_instruction',
                          'weight_hash': 'not_computed_user_instruction'})
    if not (OUT / 'task_manifest.csv').exists():
        with (OUT / 'task_manifest.csv').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(tasks[0]))
            writer.writeheader()
            writer.writerows(tasks)
    print('Prepared 30 unique Full181 tasks and two isolated YAML inputs')


def a3m_records(text):
    records = []
    for raw in text.splitlines():
        line = raw.strip().replace('\x00', '')
        if not line:
            continue
        if line.startswith('>'):
            records.append('')
        elif records:
            records[-1] += line
        else:
            raise ValueError('A3M sequence before header')
    return records


def aligned(sequence):
    return ''.join(c for c in sequence if not c.islower() and c != '.')


def msa_worker(source):
    seqs = input_sequences()
    msa = OUT / 'msa'
    raw = msa / 'colabfold_raw'
    raw.mkdir(parents=True, exist_ok=True)
    start = now()
    source = Path(source).resolve()
    spec = importlib.util.spec_from_file_location('frozen_boltz_mmseqs2', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_requests = module.requests
    count = 0

    def request(method, url, **kwargs):
        nonlocal count
        count += 1
        if count > 180:
            raise RuntimeError('MSA request ceiling reached; no unlimited service retries')
        event = {'time_utc': now(), 'request_number': count, 'method': method, 'url': url}
        try:
            response = getattr(original_requests, method)(url, **kwargs)
            event['http_status'] = response.status_code
            if 'download' not in url:
                event['response'] = response.text[:4000]
            return response
        except Exception as exc:
            event['error'] = repr(exc)
            raise
        finally:
            with (raw / 'requests.jsonl').open('a', encoding='utf-8') as handle:
                handle.write(json.dumps(event) + '\n')

    module.requests = SimpleNamespace(post=lambda url, **kw: request('post', url, **kw),
                                      get=lambda url, **kw: request('get', url, **kw))
    write_json(raw / 'batch_query.json', {'construct': 'Full181', 'condition': 'formal_msa',
               'seed': 'all_1_15', 'source_path': str(source), 'sequences': seqs,
               'pairing_scope': 'single_Full181_A_H_L_batch', 'start_utc': start,
               'boltz_source_version': '2.2.1', 'host_url': 'https://api.colabfold.com',
               'use_env': True, 'use_pairing': [True, False], 'pairing_strategy': 'greedy'})
    paired = module.run_mmseqs2(list(seqs.values()), prefix=str(raw / 'paired'),
                              use_env=True, use_pairing=True, pairing_strategy='greedy')
    unpaired = module.run_mmseqs2(list(seqs.values()), prefix=str(raw / 'unpaired'),
                                use_env=True, use_pairing=False, pairing_strategy='greedy')
    if len(paired) != 3 or len(unpaired) != 3:
        raise ValueError('expected three paired and three unpaired A3M outputs')
    pairs = [a3m_records(x) for x in paired]
    singles = [a3m_records(x) for x in unpaired]
    if len({len(p) for p in pairs}) != 1:
        raise ValueError('paired batch row counts differ; pairing keys cannot be trusted')
    manifest = []
    for i, (chain, query) in enumerate(seqs.items()):
        (raw / f'{chain}_paired.a3m').write_text(paired[i], encoding='utf-8')
        (raw / f'{chain}_unpaired.a3m').write_text(unpaired[i], encoding='utf-8')
        for family, rows in [('paired', pairs[i]), ('unpaired', singles[i])]:
            if not rows or aligned(rows[0]).replace('-', '') != query:
                raise ValueError(f'{chain} {family} query mismatch')
            if any(len(aligned(row)) != len(query) for row in rows):
                raise ValueError(f'{chain} {family} A3M aligned column mismatch')
        # Same truncation and paired/unpaired combination as the frozen generator.
        kept = [(key, sequence) for key, sequence in enumerate(pairs[i][:8192])
                if aligned(sequence) != '-' * len(query)]
        unpaired_kept = singles[i][:16384 - len(kept)]
        if kept:
            unpaired_kept = unpaired_kept[1:]
        rows = kept + [(-1, sequence) for sequence in unpaired_kept]
        path = msa / f'{chain}.csv'
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, lineterminator='\n')
            writer.writerow(['key', 'sequence'])
            writer.writerows(rows)
        manifest.append({'construct': 'Full181', 'condition': 'formal_msa', 'seed': 'all_1_15',
                         'source_path': str(path), 'chain': chain, 'query': query,
                         'query_matches': True, 'paired_raw_depth': len(pairs[i]),
                         'unpaired_raw_depth': len(singles[i]), 'combined_csv_depth': len(rows),
                         'paired_keys_retained': [key for key, _ in kept],
                         'boltz_cli_default_max_msa_seqs': 8192,
                         'actual_model_msa_depth': None,
                         'actual_model_msa_depth_note': 'unknown_until_prediction_processed_inputs_audited',
                         'pairing_scope': 'single_Full181_A_H_L_batch',
                         'fingerprint': metadata(path)})
    write_json(msa / 'msa_manifest.json', {'status': 'complete', 'started_utc': start,
               'ended_utc': now(), 'source': metadata(source), 'boltz_source_version': '2.2.1',
               'database_snapshot': 'remote_service_not_reported',
               'comparability_note': 'new search time; old database snapshot equivalence not established',
               'sha256_status': 'not_computed_user_instruction', 'chains': manifest})
    print('Complete Full181 A/H/L grouped formal MSA')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa', action='store_true')
    parser.add_argument('--msa-worker', action='store_true')
    parser.add_argument('--msa-source', type=Path, default=OUT / 'runtime_sources/boltz_2_2_1/boltz/data/msa/mmseqs2.py')
    args = parser.parse_args()
    if args.msa_worker:
        msa_worker(args.msa_source)
    elif args.msa:
        if (OUT / 'msa/msa_manifest.json').exists():
            raise RuntimeError('Existing formal MSA manifest preserved; do not silently repeat the search')
        LOG.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(Path(__file__).resolve()), '--msa-worker', '--msa-source', str(args.msa_source)]
        started = now()
        status = {'started_utc': started, 'command': command, 'timeout_seconds': 1800}
        try:
            with (LOG / 'msa_generation.stdout.log').open('w', encoding='utf-8') as stdout, (LOG / 'msa_generation.stderr.log').open('w', encoding='utf-8') as stderr:
                result = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=1800)
            status.update(exit_code=result.returncode, status='complete' if result.returncode == 0 else 'failed')
        except Exception as exc:
            status.update(exit_code=None, status='failed', error=repr(exc))
        status['ended_utc'] = now()
        write_json(OUT / 'msa_execution.json', status)
        print(json.dumps(status))
        sys.exit(0 if status['status'] == 'complete' else 1)
    else:
        prepare()


if __name__ == '__main__':
    main()
