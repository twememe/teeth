#!/usr/bin/env python3
"""Fixed Full181 formal task runner. Never runs old constructs or pilot models.

Run using the original Python environment after its resources become readable.
Artifact validation is part of accepting real formal output, not a smoke test.
"""
from __future__ import annotations
import argparse
import csv
from contextlib import contextmanager
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

from full181_contract import read_constructs, read_antibodies, build_structure_mapping
from prepare_full181_supplement import metadata, write_json, input_sequences, a3m_records, aligned

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/full181_supplement'
LOG = ROOT / 'logs/full181_supplement'


def now():
    return datetime.now(timezone.utc).isoformat()


def target_name(condition):
    return {'empty': 'full_blind', 'formal_msa': 'full_blind_msa'}[condition]


def task_dir(condition, seed):
    if condition not in ('empty', 'formal_msa') or type(seed) is not int or not 1 <= seed <= 15:
        raise ValueError('Only Full181 empty/formal_msa seeds 1-15 are authorized')
    return OUT / 'raw' / condition / f'seed{seed}'


def boltz_command(condition, seed, executable, checkpoint, cache):
    return [str(executable), 'predict', str(OUT / 'inputs' / f'{target_name(condition)}.yaml'),
            '--model', 'boltz1', '--checkpoint', str(checkpoint), '--cache', str(cache),
            '--out_dir', str(task_dir(condition, seed)), '--accelerator', 'gpu', '--devices', '1',
            '--diffusion_samples', '1', '--recycling_steps', '3', '--sampling_steps', '200',
            '--no_kernels', '--num_workers', '0', '--output_format', 'pdb',
            '--write_full_pae', '--seed', str(seed)]


def validated_input_text(condition):
    sequences = {'A': read_constructs(ROOT)['Full181'], **read_antibodies(ROOT)}
    expected = 'version: 1\nsequences:\n'
    for chain in 'AHL':
        msa = 'empty' if condition == 'empty' else f'results/full181_supplement/msa/{chain}.csv'
        expected += (f'  - protein:\n      id: {chain}\n      sequence: {sequences[chain]}\n'
                     f'      msa: {msa}\n')
    actual = (OUT / 'inputs' / f'{target_name(condition)}.yaml').read_text(encoding='utf-8')
    if actual != expected:
        raise ValueError('Full181 YAML differs from the frozen A/H/L sequence and MSA-condition contract')
    return actual


def validate_formal_msa():
    """Replay retained A3M into CSV rows; acquisition-machine metadata may differ."""
    directory = OUT / 'msa'
    raw = directory / 'colabfold_raw'
    manifest = json.loads((directory / 'msa_manifest.json').read_text(encoding='utf-8'))
    batch = json.loads((raw / 'batch_query.json').read_text(encoding='utf-8'))
    sequences = {'A': read_constructs(ROOT)['Full181'], **read_antibodies(ROOT)}
    scope = 'single_Full181_A_H_L_batch'
    if (manifest.get('status') != 'complete' or len(manifest.get('chains', [])) != 3
            or {row['chain'] for row in manifest['chains']} != set('AHL')):
        raise ValueError('Full A/H/L grouped MSA manifest is not complete')
    if (batch.get('construct') != 'Full181' or batch.get('condition') != 'formal_msa'
            or batch.get('sequences') != sequences or batch.get('pairing_scope') != scope
            or batch.get('boltz_source_version') != '2.2.1' or batch.get('use_env') is not True
            or batch.get('use_pairing') != [True, False] or batch.get('pairing_strategy') != 'greedy'):
        raise ValueError('Retained MSA batch query does not match the frozen Full181 A/H/L search')
    paired = {chain: a3m_records((raw / f'{chain}_paired.a3m').read_text(encoding='utf-8')) for chain in 'AHL'}
    unpaired = {chain: a3m_records((raw / f'{chain}_unpaired.a3m').read_text(encoding='utf-8')) for chain in 'AHL'}
    if len({len(rows) for rows in paired.values()}) != 1:
        raise ValueError('Retained A/H/L paired A3M row counts differ')
    original = {row['chain']: row for row in manifest['chains']}
    verified = []
    for chain in 'AHL':
        query = sequences[chain]
        for family, rows in (('paired', paired[chain]), ('unpaired', unpaired[chain])):
            if not rows or aligned(rows[0]).replace('-', '') != query:
                raise ValueError(f'{chain} retained {family} query mismatch')
            if any(len(aligned(row)) != len(query) for row in rows):
                raise ValueError(f'{chain} retained {family} aligned-column mismatch')
        kept = [(key, sequence) for key, sequence in enumerate(paired[chain][:8192])
                if aligned(sequence) != '-' * len(query)]
        singles = unpaired[chain][:16384-len(kept)]
        if kept:
            singles = singles[1:]
        replay = [(str(key), sequence) for key, sequence in kept] + [('-1', sequence) for sequence in singles]
        csv_path = directory / f'{chain}.csv'
        with csv_path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ['key', 'sequence']:
                raise ValueError(f'{chain} MSA CSV columns differ from the frozen format')
            actual = [(row['key'], row['sequence']) for row in reader]
        if actual != replay:
            first = next((i+1 for i, pair in enumerate(zip(actual, replay)) if pair[0] != pair[1]),
                         min(len(actual), len(replay))+1)
            raise ValueError(f'{chain} MSA CSV differs from retained A3M replay at row {first}; '
                             f'actual/replayed depths {len(actual)}/{len(replay)}')
        row = original[chain]
        if (row.get('query_matches') is not True or row.get('query') != query
                or row.get('pairing_scope') != scope or row.get('paired_raw_depth') != len(paired[chain])
                or row.get('unpaired_raw_depth') != len(unpaired[chain])
                or row.get('combined_csv_depth') != len(replay)
                or row.get('paired_keys_retained') != [key for key, _ in kept]):
            raise ValueError(f'{chain} retained MSA manifest differs from direct replay')
        verified.append({'chain': chain, 'query_matches': True, 'csv_equals_a3m_replay': True,
                         'paired_raw_depth': len(paired[chain]), 'unpaired_raw_depth': len(unpaired[chain]),
                         'combined_csv_depth': len(replay), 'current_csv_metadata': metadata(csv_path),
                         'current_paired_metadata': metadata(raw / f'{chain}_paired.a3m'),
                         'current_unpaired_metadata': metadata(raw / f'{chain}_unpaired.a3m')})
    result = {'status': 'passed', 'pairing_scope': scope, 'chains': verified,
              'current_batch_metadata': metadata(raw / 'batch_query.json'),
              'comparison_method': 'exact CSV key/sequence rows replayed from retained A3M; no hashes',
              'acquisition_metadata_used_as_deployment_identity': False}
    write_json(directory / 'current_input_validation.json', {'validated_utc': now(), **result})
    return result


def validate_saved_state(condition, seed, state, artifacts):
    if (state.get('construct'), state.get('condition'), state.get('seed')) != ('Full181', condition, seed):
        raise ValueError('Per-task state identity does not match the requested condition/seed')
    provenance = state.get('provenance', {})
    if provenance.get('input_text') != validated_input_text(condition):
        raise ValueError('Saved task input differs from the current validated input')
    if provenance.get('sequences') != {'A': read_constructs(ROOT)['Full181'], **read_antibodies(ROOT)}:
        raise ValueError('Saved task antigen/antibody sequence identity mismatch')
    for key in ('checkpoint', 'ccd'):
        saved = provenance.get(key, {})
        if not saved.get('exists') or metadata(saved['path']) != saved:
            raise ValueError(f'Saved task {key} metadata does not match its current source')
    executable = provenance.get('boltz_executable', {})
    if not executable.get('exists') or metadata(executable['path']) != executable:
        raise ValueError('Saved Boltz executable metadata changed')
    expected_command = boltz_command(condition, seed, executable['path'],
                                     provenance['checkpoint']['path'], provenance.get('cache_path', ''))
    if provenance.get('command') != expected_command:
        raise ValueError('Saved command does not match this condition/seed and the frozen formal parameters')
    saved_msa = provenance.get('msa', [])
    if len(saved_msa) != (3 if condition == 'formal_msa' else 0):
        raise ValueError('Saved MSA condition provenance mismatch')
    for saved in saved_msa:
        if not saved.get('exists') or metadata(saved['path']) != saved:
            raise ValueError('Saved MSA source metadata changed')
    if state.get('artifact_fingerprints') != artifacts['artifact_fingerprints']:
        raise ValueError('Accepted artifact metadata changed; preserve the task for review')


def validate_artifacts(condition, seed):
    import numpy as np
    target = target_name(condition)
    directory = task_dir(condition, seed) / f'boltz_results_{target}' / 'predictions' / target
    pdb = directory / f'{target}_model_0.pdb'
    confidence = directory / f'confidence_{target}_model_0.json'
    pae = directory / f'pae_{target}_model_0.npz'
    # This validates exact Full181/VH119/VL106 identity, not only chain lengths.
    mapping = build_structure_mapping(pdb, 'Full181', ROOT, require_complete=True)
    lengths = {chain: mapping[chain]['metadata']['observed_length'] for chain in 'AHL'}
    if lengths != {'A': 181, 'H': 119, 'L': 106}:
        raise ValueError(f'Unexpected verified chain lengths: {lengths}')
    atom_count = 0
    model_count = 0
    b_factors = []
    for line in pdb.read_text(encoding='utf-8').splitlines():
        if line.startswith('MODEL '):
            model_count += 1
        if line.startswith('ATOM  '):
            atom_count += 1
            xyz = [float(line[start:start+8]) for start in (30, 38, 46)]
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError('PDB contains nonfinite ATOM coordinates')
            value = float(line[60:66])
            if not math.isfinite(value) or not 0 <= value <= 100:
                raise ValueError('PDB pLDDT field is not finite on the frozen 0-100 scale')
            b_factors.append(value)
    if not atom_count or model_count > 1:
        raise ValueError('Exactly one PDB coordinate model is required')
    metrics = json.loads(confidence.read_text(encoding='utf-8'))
    for key in ('confidence_score', 'ptm', 'iptm', 'complex_plddt'):
        value = metrics.get(key)
        if type(value) not in (float, int) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'Missing, nonfinite, or out-of-range 0-1 confidence field: {key}')
    with np.load(pae, allow_pickle=False) as data:
        matrix = data['pae']
        if (matrix.shape != (406, 406) or matrix.dtype.kind not in 'fiu'
                or not np.isfinite(matrix).all() or (matrix < 0).any()):
            raise ValueError(f'Invalid PAE array: {matrix.shape}')
        pae_metadata = {'shape': list(matrix.shape), 'dtype': str(matrix.dtype),
                        'finite': True, 'minimum': float(matrix.min()),
                        'maximum': float(matrix.max()), 'units': 'angstrom'}
    return {'pdb_path': str(pdb), 'confidence_path': str(confidence), 'pae_path': str(pae),
            'chain_lengths': lengths, 'validation': 'passed', 'validated_utc': now(),
            'sequence_validation': {c: mapping[c]['metadata'] for c in 'AHL'},
            'pae_validation': pae_metadata, 'confidence_json_scale': '0-1',
            'pdb_plddt_scale': '0-100 per frozen Boltz PDB writer',
            'pdb_plddt_observed_range': [min(b_factors), max(b_factors)],
            'artifact_fingerprints': [metadata(path) for path in (pdb, confidence, pae)]}


@contextmanager
def manifest_mutex():
    """Serialize shared manifest snapshots; the OS releases locks on process exit."""
    lock_path = OUT / 'task_manifest.lock'
    with lock_path.open('a+b') as handle:
        if os.name == 'nt':
            import msvcrt
            if handle.tell() == 0:
                handle.write(b'0')
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def refresh_manifest(block_reason=None):
    with manifest_mutex():
        return _refresh_manifest(block_reason)


def _refresh_manifest(block_reason=None):
    rows = []
    for condition in ('empty', 'formal_msa'):
        for seed in range(1, 16):
            directory = task_dir(condition, seed)
            state_path = directory / 'run_state.json'
            state = json.loads(state_path.read_text(encoding='utf-8')) if state_path.exists() else {}
            rows.append({'construct': 'Full181', 'condition': condition, 'seed': seed,
                         'task_id': f'Full181_{condition}_seed{seed}', 'source_path': str(directory),
                         'input_path': str(OUT / 'inputs' / f'{target_name(condition)}.yaml'),
                         'status': state.get('status', 'blocked_resource' if block_reason else 'planned'),
                         'failure_reason': state.get('failure_reason', block_reason or ''),
                         'exit_code': state.get('exit_code', ''),
                         'started_utc': state.get('started_utc', ''), 'ended_utc': state.get('ended_utc', ''),
                         'input_hash': 'not_computed_user_instruction',
                         'weight_hash': 'not_computed_user_instruction',
                         'provenance_path': str(state_path) if state_path.exists() else ''})
    path = OUT / 'task_manifest.csv'
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='', dir=path.parent,
                                     prefix=path.name + '.tmp_', delete=False) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)
    return rows


def environment_audit(boltz, checkpoint, cache):
    record = {'time_utc': now(), 'python': sys.version, 'executable': sys.executable,
              'checkpoint': metadata(checkpoint), 'ccd': metadata(cache / 'ccd.pkl'),
              'boltz_executable': metadata(boltz), 'hash_status': 'not_computed_user_instruction'}
    for name in ('boltz', 'torch', 'numpy', 'biopython', 'prodigy-prot', 'pytorch-lightning'):
        try:
            record[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            record[name] = None
    write_json(OUT / 'environment_fingerprint.json', record)
    if record['boltz'] != '2.2.1' or record['torch'] != '2.6.0+cu124':
        raise RuntimeError('Frozen Boltz 2.2.1 / torch 2.6.0+cu124 environment unavailable')
    if boltz.parent.resolve() != Path(sys.executable).parent.resolve():
        raise RuntimeError('Boltz executable must belong to the currently audited original environment')
    if not checkpoint.is_file() or checkpoint.name != 'boltz1_conf.ckpt':
        raise RuntimeError('Original boltz1_conf.ckpt is unavailable')
    if checkpoint.stat().st_size != 3595352714:
        raise RuntimeError('Checkpoint size differs from historical inventory; no model substitution')
    if not (cache / 'ccd.pkl').is_file() or (cache / 'ccd.pkl').stat().st_size != 345859128:
        raise RuntimeError('Original CCD cache missing or differs in size from historical inventory')
    help_result = subprocess.run([str(boltz), 'predict', '--help'], capture_output=True, text=True)
    (LOG / 'boltz_predict_help.txt').write_text(help_result.stdout + help_result.stderr, encoding='utf-8')
    if help_result.returncode:
        raise RuntimeError('Frozen environment CLI cannot be read')
    for flag in ('--model', '--checkpoint', '--cache', '--seed', '--out_dir', '--accelerator', '--devices',
                 '--diffusion_samples', '--recycling_steps', '--sampling_steps', '--no_kernels',
                 '--num_workers', '--output_format', '--write_full_pae'):
        if flag not in help_result.stdout:
            raise RuntimeError(f'Frozen CLI missing required flag: {flag}')
    return record


def gpu_snapshot(gpu):
    result = subprocess.run(['nvidia-smi', '-i', str(gpu),
                             '--query-gpu=index,name,memory.total,memory.used,utilization.gpu',
                             '--format=csv,noheader,nounits'], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip())
    fields = [x.strip() for x in result.stdout.strip().split(',')]
    return {'time_utc': now(), 'gpu': fields[0], 'name': fields[1],
            'memory_total_mib': int(fields[2]), 'memory_used_mib': int(fields[3]),
            'utilization_percent': int(fields[4])}


def run_branch(args):
    args.boltz, args.checkpoint, args.cache = (path.resolve() for path in (args.boltz, args.checkpoint, args.cache))
    read_constructs(ROOT)
    read_antibodies(ROOT)
    input_sequences()
    LOG.mkdir(parents=True, exist_ok=True)
    environment = environment_audit(args.boltz, args.checkpoint, args.cache)
    validated_input_text(args.condition)
    prodigy = args.prodigy or Path(sys.executable).parent / 'prodigy'
    analysis_script = ROOT / 'src/analyze_full181_supplement.py'
    if not analysis_script.is_file() or not prodigy.is_file() or environment['prodigy-prot'] != '2.4.0':
        raise RuntimeError('The original PRODIGY 2.4.0 and contact-analysis entry point are required before seed 1')
    if prodigy.parent.resolve() != Path(sys.executable).parent.resolve():
        raise RuntimeError('PRODIGY executable must belong to the audited original environment')
    msa_validation = validate_formal_msa() if args.condition == 'formal_msa' else None
    for seed in range(1, 16):
        directory = task_dir(args.condition, seed)
        directory.mkdir(parents=True, exist_ok=True)
        yaml = OUT / 'inputs' / f'{target_name(args.condition)}.yaml'
        command = boltz_command(args.condition, seed, args.boltz, args.checkpoint, args.cache)
        provenance = {'command': command, 'input_text': validated_input_text(args.condition),
                      'input_file': metadata(yaml),
                      'boltz_executable': metadata(args.boltz),
                      'cache_path': str(args.cache),
                      'sequences': input_sequences(), 'checkpoint': metadata(args.checkpoint),
                      'ccd': metadata(args.cache / 'ccd.pkl'),
                      'msa': [metadata(OUT / 'msa' / f'{c}.csv') for c in 'AHL'] if args.condition == 'formal_msa' else [],
                      'msa_replay': msa_validation,
                      'environment': {key: environment[key] for key in ('boltz', 'torch', 'python', 'executable')}}
        state_path = directory / 'run_state.json'
        if state_path.exists():
            previous = json.loads(state_path.read_text(encoding='utf-8'))
            if previous.get('status') in ('complete', 'downstream_blocked') and previous.get('provenance') == provenance:
                artifacts = validate_artifacts(args.condition, seed)
                validate_saved_state(args.condition, seed, previous, artifacts)
                if seed == 1:
                    previous.update(status='complete', downstream_status='pending')
                    write_json(state_path, previous)
                    try:
                        previous['downstream'] = run_seed1_downstream(args.condition, prodigy, artifacts)
                        previous['downstream_status'] = 'complete'
                    except Exception as exc:
                        previous.update(status='downstream_blocked', downstream_status='blocked', failure_reason=repr(exc))
                        write_json(state_path, previous)
                        refresh_manifest()
                        raise
                    write_json(state_path, previous)
                print(f'Validated existing {args.condition} seed {seed}', flush=True)
                continue
            raise RuntimeError(f'Existing incomplete/mismatched task preserved: {directory}; no automatic rerun')
        snapshot = gpu_snapshot(args.gpu)
        if snapshot['memory_used_mib'] > 1024 or snapshot['utilization_percent'] > 5:
            raise RuntimeError(f'GPU {args.gpu} is occupied; no processes terminated')
        free = shutil.disk_usage(directory).free
        if free < 2 * 1024**3:
            raise RuntimeError(f'Insufficient task filesystem free space: {free} bytes')
        state = {'construct': 'Full181', 'condition': args.condition, 'seed': seed,
                 'source_path': str(directory), 'status': 'running', 'started_utc': now(),
                 'gpu': args.gpu, 'provenance': provenance}
        write_json(state_path, state)
        log_dir = LOG / args.condition / f'seed{seed}'
        log_dir.mkdir(parents=True, exist_ok=True)
        write_json(log_dir / 'command.json', command)
        env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(args.gpu))
        started = time.monotonic()
        try:
            with (log_dir / 'stdout.log').open('w', encoding='utf-8') as stdout, (log_dir / 'stderr.log').open('w', encoding='utf-8') as stderr, (log_dir / 'gpu.jsonl').open('w', encoding='utf-8') as gpu_log:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
                while process.poll() is None:
                    try:
                        sample = gpu_snapshot(args.gpu)
                    except Exception as exc:
                        # A monitoring failure must not orphan an active formal prediction.
                        sample = {'time_utc': now(), 'monitor_error': repr(exc)}
                    gpu_log.write(json.dumps(sample) + '\n')
                    gpu_log.flush()
                    time.sleep(10)
            state['exit_code'] = process.returncode
            if process.returncode:
                raise RuntimeError(f'Boltz exited {process.returncode}; see {log_dir}')
            artifacts = validate_artifacts(args.condition, seed)
            state.update(artifacts, status='complete', prediction_status='complete', failure_reason='')
            if seed == 1:
                state['downstream_status'] = 'pending'
                write_json(state_path, state)
                state['downstream'] = run_seed1_downstream(args.condition, prodigy, artifacts)
                state['downstream_status'] = 'complete'
        except Exception as exc:
            state.update(status='downstream_blocked' if state.get('prediction_status') == 'complete' else 'failed',
                         failure_reason=repr(exc))
        state.update(ended_utc=now(), wall_seconds=time.monotonic() - started)
        write_json(state_path, state)
        refresh_manifest()
        print(f'{args.condition} seed {seed}: {state["status"]}', flush=True)
        if state['status'] != 'complete':
            # Seed 1 must complete the formal artifact chain before continuing.
            if seed == 1:
                raise RuntimeError(state['failure_reason'])
    refresh_manifest()


def prodigy_row(condition, seed, source_path=''):
    return {'construct': 'Full181', 'condition': condition, 'seed': seed,
            'source_path': source_path or str(task_dir(condition, seed)), 'status': 'blocked_prediction',
            'failure_reason': 'No accepted Full181 prediction', 'exit_code': '',
            'prodigy_dG_kcal_mol': '', 'prodigy_Kd_M_25C': '',
            'temperature_C': '', 'configured_temperature_C': 25, 'intermolecular_contacts': '',
            'stdout_path': '', 'stderr_path': '', 'command_path': '',
            'started_utc': '', 'ended_utc': '',
            'method_note': 'computational proxy; original unrelaxed A versus H,L; not measured KD'}


def run_prodigy_one(executable, condition, seed, artifacts):
    version = importlib.metadata.version('prodigy-prot')
    if version != '2.4.0':
        raise RuntimeError(f'Expected original PRODIGY 2.4.0, found {version}')
    row = prodigy_row(condition, seed, artifacts['pdb_path'])
    log = LOG / condition / f'seed{seed}' / 'prodigy'
    log.mkdir(parents=True, exist_ok=True)
    command = [str(executable), artifacts['pdb_path'], '--selection', 'A', 'H,L', '--temperature', '25']
    record = {'construct': 'Full181', 'condition': condition, 'seed': seed,
              'source_path': artifacts['pdb_path'], 'command': command, 'version': version,
              'artifact_fingerprints': artifacts['artifact_fingerprints']}
    existing_path = task_dir(condition, seed) / 'prodigy_result.json'
    if existing_path.is_file():
        existing = json.loads(existing_path.read_text(encoding='utf-8'))
        command_path = log / 'command.json'
        saved_command = json.loads(command_path.read_text(encoding='utf-8')) if command_path.is_file() else None
        if saved_command != record:
            raise RuntimeError('Existing PRODIGY provenance differs; original attempt retained')
        if not all((log / name).is_file() for name in ('stdout.log', 'stderr.log')):
            raise RuntimeError('Existing PRODIGY raw log is missing; original attempt retained')
        return existing
    write_json(log / 'command.json', record)
    row.update(started_utc=now(), stdout_path=str(log / 'stdout.log'),
               stderr_path=str(log / 'stderr.log'), command_path=str(log / 'command.json'))
    try:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        (log / 'stdout.log').write_text(result.stdout, encoding='utf-8')
        (log / 'stderr.log').write_text(result.stderr, encoding='utf-8')
        row['exit_code'] = result.returncode
        dg = re.search(r'binding affinity .*?:\s+([-+0-9.eE]+)', result.stdout, re.I)
        kd = re.search(r'dissociation constant .*?:\s+([-+0-9.eE]+)', result.stdout, re.I)
        nc = re.search(r'No\. of intermolecular contacts:\s+(\d+)', result.stdout)
        temperature = re.search(r'dissociation constant \(M\) at\s+([-+0-9.]+)', result.stdout, re.I)
        # The exact no-contact exception is from the frozen PRODIGY 2.4.0 source.
        # It is an observed absence of an interface, not a reason to reject a seed.
        if 'No contacts found for selection' in result.stdout + result.stderr or (nc and int(nc[1]) == 0):
            row.update(status='no_interface', failure_reason='No contacts for A versus H,L at PRODIGY cutoff',
                       intermolecular_contacts=0)
        elif result.returncode == 0 and dg and kd and nc and temperature:
            delta_g, dissociation, actual_temperature = float(dg[1]), float(kd[1]), float(temperature[1])
            if not all(math.isfinite(x) for x in (delta_g, dissociation, actual_temperature)):
                raise ValueError('Nonfinite PRODIGY output')
            if dissociation <= 0 or actual_temperature != 25 or int(nc[1]) <= 0:
                raise ValueError('Invalid PRODIGY Kd, temperature, or contact count')
            row.update(status='complete', failure_reason='', prodigy_dG_kcal_mol=delta_g,
                       prodigy_Kd_M_25C=dissociation, intermolecular_contacts=int(nc[1]),
                       temperature_C=actual_temperature)
        else:
            row.update(status='failed', failure_reason='Nonzero exit or unparseable required PRODIGY fields; see logs')
    except Exception as exc:
        row.update(status='failed', failure_reason=repr(exc))
        if not (log / 'stderr.log').exists():
            (log / 'stderr.log').write_text(repr(exc) + '\n', encoding='utf-8')
    row['ended_utc'] = now()
    write_json(task_dir(condition, seed) / 'prodigy_result.json', row)
    return row


def run_seed1_downstream(condition, executable, artifacts):
    """Complete real formal seed-1 analysis before advancing to seeds 2-15."""
    refresh_manifest()
    log = LOG / condition / 'seed1' / 'analysis'
    log.mkdir(parents=True, exist_ok=True)
    script = ROOT / 'src/analyze_full181_supplement.py'
    analysis_out = task_dir(condition, 1) / 'analysis'
    command = [sys.executable, '-B', str(script), '--condition', condition, '--seed', '1', '--output', str(analysis_out)]
    write_json(log / 'command.json', {'command': command, 'script': metadata(script),
               'construct': 'Full181', 'condition': condition, 'seed': 1,
               'source_path': artifacts['pdb_path']})
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    (log / 'stdout.log').write_text(result.stdout, encoding='utf-8')
    (log / 'stderr.log').write_text(result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'Formal seed-1 contact analysis failed: {log}')
    summary_path, contact_path = analysis_out / 'full181_model_summary.csv', analysis_out / 'full181_interface_contacts.csv'
    with summary_path.open(encoding='utf-8-sig', newline='') as handle:
        summaries = [row for row in csv.DictReader(handle)
                     if row['construct'] == 'Full181' and row['condition'] == condition and int(row['seed']) == 1]
    if len(summaries) != 1:
        raise RuntimeError('Formal seed-1 contact analysis did not return exactly one model summary')
    source = Path(summaries[0]['source_path'])
    source = source if source.is_absolute() else ROOT / source
    if source.resolve() != Path(artifacts['pdb_path']).resolve():
        raise RuntimeError('Formal seed-1 contact summary refers to a different source PDB')
    with contact_path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        if not {'construct', 'condition', 'seed', 'source_path'}.issubset(reader.fieldnames or []):
            raise RuntimeError('Contact table does not carry per-model provenance')
        count = 0
        for row in reader:
            if row['construct'] != 'Full181' or row['condition'] != condition or int(row['seed']) != 1:
                continue
            contact_source = Path(row['source_path'])
            contact_source = contact_source if contact_source.is_absolute() else ROOT / contact_source
            if contact_source.resolve() != Path(artifacts['pdb_path']).resolve():
                raise RuntimeError('Formal seed-1 contacts refer to a different source PDB')
            count += 1
    score = run_prodigy_one(executable, condition, 1, artifacts)
    if score['status'] not in ('complete', 'no_interface'):
        raise RuntimeError('Formal seed-1 PRODIGY execution/parsing failed; raw model is preserved')
    return {'status': 'complete', 'source_path': artifacts['pdb_path'],
            'contact_summary_path': str(summary_path), 'contact_summary': summaries[0],
            'contact_table_path': str(contact_path), 'contact_rows_at_5A': count,
            'prodigy_status': score['status'], 'prodigy_result_path': str(task_dir(condition, 1) / 'prodigy_result.json'),
            'no_contact_models_are_retained': True, 'completed_utc': now()}


def run_prodigy(executable):
    # This runs only on accepted original, unrelaxed Full formal outputs.
    rows = []
    for condition in ('empty', 'formal_msa'):
        for seed in range(1, 16):
            state_path = task_dir(condition, seed) / 'run_state.json'
            state = json.loads(state_path.read_text(encoding='utf-8')) if state_path.exists() else {}
            row = prodigy_row(condition, seed, state.get('pdb_path', ''))
            if state.get('status') in ('complete', 'downstream_blocked'):
                try:
                    artifacts = validate_artifacts(condition, seed)
                    validate_saved_state(condition, seed, state, artifacts)
                    row = run_prodigy_one(executable, condition, seed, artifacts)
                except Exception as exc:
                    row.update(status='failed_validation', failure_reason=repr(exc))
            rows.append(row)
    with (OUT / 'full181_prodigy.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--condition', choices=['empty', 'formal_msa'])
    parser.add_argument('--gpu', type=int)
    parser.add_argument('--boltz', type=Path)
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--refresh-manifest', action='store_true')
    parser.add_argument('--block-reason')
    parser.add_argument('--prodigy', type=Path)
    args = parser.parse_args()
    if args.refresh_manifest:
        refresh_manifest(args.block_reason)
    elif all(x is not None for x in (args.condition, args.gpu, args.boltz, args.checkpoint, args.cache)):
        run_branch(args)
    elif args.prodigy:
        run_prodigy(args.prodigy)
    else:
        parser.error('Provide condition, gpu, boltz, checkpoint, and cache; or a reporting mode')


if __name__ == '__main__':
    main()
