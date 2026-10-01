#!/usr/bin/env python3
"""Export current supplement provenance and schemas without hashes or inference."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import importlib.metadata

from prepare_full181_supplement import metadata, write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/full181_supplement'
REPORT = ROOT / 'reports/full181_supplement'
LOG = ROOT / 'logs/full181_supplement'


def csv_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f'No rows for {path}')
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    timestamp = datetime.now(timezone.utc).isoformat()
    tasks = csv_rows(OUT / 'task_manifest.csv')
    unique = {(r['construct'], r['condition'], int(r['seed'])) for r in tasks}
    expected = {('Full181', c, s) for c in ('empty', 'formal_msa') for s in range(1, 16)}
    if len(tasks) != 30 or unique != expected:
        raise ValueError('Task manifest does not contain exactly the authorized 30 unique tasks')
    seq = json.loads((OUT / 'sequence_audit/sequence_audit.json').read_text(encoding='utf-8'))
    msa = json.loads((OUT / 'msa/msa_manifest.json').read_text(encoding='utf-8'))
    resource = json.loads((OUT / 'resource_audit.json').read_text(encoding='utf-8'))
    completed = sum(r['status'] == 'complete' for r in tasks)
    summaries = csv_rows(OUT / 'full181_model_summary.csv')
    prodigy = csv_rows(OUT / 'full181_prodigy.csv')
    analyses = sum(r.get('status') == 'complete' for r in summaries)
    affinity_complete = sum(r.get('status') == 'complete' for r in prodigy)
    no_interface = sum(r.get('status') == 'no_interface' for r in prodigy)
    local_raw_counts = {label: len(list((OUT / 'raw').rglob(pattern))) for label, pattern in
                        [('pdb', '*_model_0.pdb'), ('confidence', 'confidence_*.json'), ('pae', 'pae_*.npz')]}
    for label, rows in [('summary', summaries), ('prodigy', prodigy)]:
        keys = {(r['construct'], r['condition'], int(r['seed'])) for r in rows}
        if len(rows) != 30 or keys != expected:
            raise ValueError(f'{label} must retain exactly 30 authorized keys')
    write_json(OUT / 'input_audit.json', {'time_utc': timestamp, 'sequence_audit': seq,
               'resource_audit': resource, 'msa_status': msa['status'],
               'hash_policy': 'not_computed_per_direct_user_instruction',
               'formal_models_completed': completed,
               'formal_models_planned': 30})
    inputs = [ROOT / 'data' / n for n in ('IL24_full_1_181.fasta', 'IL24_native_27_181.fasta',
               'IL24_immunogen_27_160.fasta', 'VH_VL_clean.fasta', 'antibody_imgt_numbering.csv',
               'cdr_annotation.csv', 'il24_numbering_map.csv')]
    inputs += [p for p in (OUT / 'inputs').iterdir() if p.is_file()]
    inputs += [OUT / 'msa' / f'{chain}.csv' for chain in 'AHL']
    input_rows = []
    for path in inputs:
        stat = metadata(path)
        input_rows.append({'construct': 'input_contract', 'condition': 'input_contract', 'seed': 'NA',
                           'source_path': stat['path'], 'bytes': stat.get('bytes', 'NA'),
                           'mtime_ns': stat.get('mtime_ns', 'NA'), 'hash_status': 'not_computed_user_instruction'})
    write_csv(OUT / 'input_fingerprints.csv', input_rows)
    previous = json.loads((OUT / 'baseline_file_metadata.json').read_text(encoding='utf-8'))
    changed = []
    for row in previous:
        current = metadata(row['path'])
        if any(row.get(k) != current.get(k) for k in ('exists', 'bytes', 'mtime_ns')):
            changed.append({'before': row, 'after': current})
    verification = {'time_utc': timestamp, 'task_count': len(tasks), 'unique_task_count': len(unique),
                    'completed_models': completed,
                    'completed_contact_analyses': analyses,
                    'numeric_prodigy_observations': affinity_complete,
                    'prodigy_no_interface_observations': no_interface,
                    'local_original_raw_file_counts': local_raw_counts,
                    'blocked_tasks': sum(r['status'].startswith('blocked') for r in tasks),
                    'old_files_metadata_checked': len(previous), 'old_files_metadata_changed': changed,
                    'old_file_content_hash_verified': False,
                    'hash_note': 'No hash computed; unchanged metadata is not a cryptographic content guarantee',
                    'formal_msa_status': msa['status'],
                    'formal_msa_query_checks': {r['chain']: r['query_matches'] for r in msa['chains']},
                    'mapping_test_log': str(OUT / 'sequence_audit/mapping_tests.log'),
                    'full_prediction_acceptance': 'passed_30_formal_predictions' if completed == 30 and all(n == 30 for n in local_raw_counts.values()) else 'incomplete_see_per_task_status',
                    'full_downstream_acceptance': 'passed' if analyses == 30 and affinity_complete + no_interface == 30 else 'incomplete_see_per_task_status'}
    for name in ('server_artifact_validation.json', 'server_baseline_metadata_verification.json', 'local_transfer_validation.json'):
        path = OUT / name
        if path.is_file():
            verification[name.removesuffix('.json')] = json.loads(path.read_text(encoding='utf-8'))
    write_json(OUT / 'delivery_verification.json', verification)
    runtime = {'python': sys.version, 'executable': sys.executable, 'role': 'non_GPU_preparation_and_data_export',
               'gpu_environment': resource.get('original_environment', 'see_resource_audit')}
    for package in ('numpy', 'pandas', 'requests'):
        try:
            runtime[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            runtime[package] = None
    write_json(OUT / 'preparation_environment.json', runtime)
    schema = []
    for path in sorted(OUT.rglob('*.csv')):
        if any(x in path.parts for x in ('runtime_sources', 'server_baseline_snapshot')) or path.name == 'artifact_fingerprints.csv':
            continue
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.reader(handle)
            fields = next(reader, [])
            count = sum(1 for _ in reader)
        schema.append({'path': path.relative_to(ROOT).as_posix(), 'rows': count, 'fields': fields})
    write_json(OUT / 'table_schemas.json', schema)
    artifacts = [p for folder in (OUT, LOG, REPORT) for p in folder.rglob('*')
                 if p.is_file() and p.name != 'artifact_fingerprints.csv' and '__pycache__' not in p.parts
                 and not p.name.startswith('Full181_supplement_delivery')]
    artifacts += list((ROOT / 'src').glob('*full181*.py')) + [ROOT / 'tests/test_full181_mapping.py']
    rows = []
    for path in sorted(set(artifacts)):
        record = metadata(path)
        rows.append({'source_path': record['path'], 'bytes': record['bytes'],
                     'mtime_ns': record['mtime_ns'], 'hash_status': 'not_computed_user_instruction'})
    write_csv(OUT / 'artifact_fingerprints.csv', rows)
    print(json.dumps({k: verification[k] for k in ('task_count', 'unique_task_count', 'completed_models',
                     'blocked_tasks', 'old_files_metadata_checked', 'old_files_metadata_changed', 'formal_msa_status')}, ensure_ascii=True))
    print(f'Exported schemas for {len(schema)} tables and metadata for {len(rows)} files')


if __name__ == '__main__':
    main()
