#!/usr/bin/env python3
"""Audit accepted formal artifacts and actual stored MSA arrays; no inference."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from prepare_full181_supplement import metadata, write_json
from run_full181_supplement import ROOT, OUT, LOG, now, validate_artifacts, validate_saved_state


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    from boltz.data.parse.csv import parse_csv
    expected = {c: parse_csv(OUT / 'msa' / f'{c}.csv', max_seqs=8192) for c in 'AHL'}
    acquisition = json.loads((OUT / 'msa/msa_manifest.json').read_text())
    depths = {r['chain']: r for r in acquisition['chains']}
    task_rows, msa_rows, failures = [], [], []
    for condition in ('empty', 'formal_msa'):
        for seed in range(1, 16):
            directory = OUT / 'raw' / condition / f'seed{seed}'
            state_path = directory / 'run_state.json'
            row = dict(construct='Full181', condition=condition, seed=seed, source_path=str(state_path))
            try:
                state = json.loads(state_path.read_text())
                if state.get('status') not in ('complete', 'downstream_blocked'):
                    raise ValueError(f"Task status {state.get('status')}")
                artifacts = validate_artifacts(condition, seed)
                validate_saved_state(condition, seed, state, artifacts)
                row.update(status='passed', prediction_state=state['status'], failure_reason='',
                           pdb_path=artifacts['pdb_path'], confidence_path=artifacts['confidence_path'],
                           pae_path=artifacts['pae_path'], wall_seconds=state['wall_seconds'],
                           started_utc=state['started_utc'], ended_utc=state['ended_utc'], exit_code=state['exit_code'])
                samples = [json.loads(line) for line in (LOG / condition / f'seed{seed}/gpu.jsonl').read_text().splitlines()]
                valid_samples = [s for s in samples if 'memory_used_mib' in s]
                row['gpu_index'] = state['gpu']
                row['gpu_peak_observed_memory_mib'] = max(s['memory_used_mib'] for s in valid_samples) if valid_samples else 'NA'
                row['gpu_monitor_samples'] = len(valid_samples)
                prediction = Path(artifacts['pdb_path']).parent
                processed = prediction.parents[1] / 'processed'
                record_path = processed / 'records' / f'{prediction.name}.json'
                record = json.loads(record_path.read_text())
                chains = record['chains']
                if {c['chain_name']:c['num_residues'] for c in chains} != {'A':181,'H':119,'L':106}:
                    raise ValueError('Processed record chain identity differs')
                for chain in chains:
                    name = chain['chain_name']
                    current = dict(construct='Full181', condition=condition, seed=seed,
                                   source_path=str(record_path), chain=name, cli_max_msa_seqs=8192,
                                   actual_model_tensor_depth='unknown_not_instrumented',
                                   model_tensor_depth_not_equated_to_stored_depth=True)
                    if condition == 'empty':
                        if chain['msa_id'] != -1:
                            raise ValueError('Empty condition has a nonempty stored MSA ID')
                        current.update(status='passed', stored_msa_depth='NA', msa_id=-1,
                                       configured_mode='explicit_empty_query_only', csv_replay_matches='NA')
                    else:
                        path = processed / 'msa' / f"{chain['msa_id']}.npz"
                        with np.load(path, allow_pickle=False) as saved:
                            matches = all(np.array_equal(saved[k], getattr(expected[name], k))
                                          for k in ('sequences', 'residues', 'deletions'))
                            actual_depth = len(saved['sequences'])
                        if not matches:
                            raise ValueError(f'{name} processed MSA differs from this batch CSV parser replay')
                        d = depths[name]
                        current.update(status='passed', source_path=str(path), stored_msa_depth=actual_depth,
                                       paired_raw_depth=d['paired_raw_depth'], unpaired_raw_depth=d['unpaired_raw_depth'],
                                       final_csv_depth=d['combined_csv_depth'], csv_replay_matches=True,
                                       processing_rule='original Boltz 2.2.1 parse_csv deduplication then 8192 cap',
                                       msa_id=chain['msa_id'], paired_key_origin='single_Full181_A_H_L_batch')
                    msa_rows.append(current)
            except Exception as exc:
                row.update(status='failed_or_missing', failure_reason=repr(exc))
                failures.append(row)
            task_rows.append(row)
    write_csv(OUT / 'execution_per_seed_audit.csv', task_rows)
    if msa_rows:
        write_csv(OUT / 'processed_msa_audit.csv', msa_rows)
    result = dict(audited_utc=now(), expected_tasks=30, validated_tasks=30-len(failures),
                  failed_or_missing=failures, hashes_computed=False, inference_performed=False,
                  method='existing artifact/sequence/parameter validation plus exact stored MSA array comparison')
    write_json(OUT / 'server_artifact_validation.json', result)
    previous = json.loads((OUT / 'server_baseline_file_metadata.json').read_text())
    changed = []
    for old in previous:
        current = metadata(old['path'])
        if not current.get('exists') or any(old[k] != current[k] for k in ('bytes', 'mtime_ns')):
            changed.append({'before':old,'after':current})
    write_json(OUT / 'server_baseline_metadata_verification.json',
               dict(audited_utc=now(), checked_files=len(previous), changed_files=changed,
                    hash_verified=False, method='path/size/mtime; not cryptographic equality'))
    print(json.dumps(dict(validated_tasks=30-len(failures), processed_msa_rows=len(msa_rows),
                          original_files_checked=len(previous), original_metadata_changes=len(changed))))
    if changed or args.require_complete and failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
