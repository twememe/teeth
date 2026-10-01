#!/usr/bin/env python3
"""Audit the 60 snapshotted legacy formal models on CPU; never run predictions.

Only Native155/Immunogen134, empty/formal_msa, and seeds 1-15 are addressed.
Source snapshots are read-only. Outputs are flushed after each model, with a
JSON checkpoint after every five models. No digests are computed.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import traceback

import numpy as np

from analyze_full181_supplement import analyze_model, SUMMARY, CONTACTS, GEOMETRY, CDR
from full181_contract import build_structure_mapping, read_antibodies, read_constructs

ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTS = {'Native155': ('native_blind', 155), 'Immunogen134': ('immunogen_blind', 134)}
CONDITIONS = {'empty': 'formal', 'formal_msa': 'formal_msa'}
SCORES = ('confidence_score', 'ptm', 'iptm', 'complex_plddt')
DISTANCE_TOLERANCE_A = 0.00055  # Saved distances have 3 decimals; legacy coordinates used float32.
EXTRA = ['remote_source_path', 'server_source_path', 'source_inventory_path', 'source_condition',
         'raw_validation_status', 'saved_comparison_status', 'evidence_level']
VALIDATION = ['construct', 'condition', 'seed', 'source_path', *EXTRA,
              'status', 'failure_reason', 'chain_A_length', 'chain_H_length', 'chain_L_length',
              'sequence_identity', 'confidence_valid', 'pae_valid', 'pae_shape', 'pae_min_A', 'pae_max_A',
              'processed_record_valid', 'chain_id_mapping', 'chain_pair_confidence_valid',
              'source_inventory_bytes_match', 'saved_summary_matches', 'summary_difference_count',
              'saved_contact_pairs', 'recomputed_contact_pairs', 'missing_recomputed_contact_pairs',
              'added_recomputed_contact_pairs', 'contact_field_difference_count',
              'distance_outside_saved_precision_count', 'distance_rounded_3dp_difference_count',
              'distance_max_abs_difference_A', 'distance_abs_tolerance_A',
              'paratope_set_matches_saved_contacts', 'cdr_counts_match_saved_contacts',
              'confidence_path', 'remote_confidence_path', 'pae_path', 'remote_pae_path',
              'processed_record_path', 'remote_processed_record_path',
              'saved_summary_path', 'remote_saved_summary_path', 'saved_contacts_path', 'remote_saved_contacts_path']


def now():
    return datetime.now(timezone.utc).isoformat()


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def residue_set(text):
    return {int(x) for x in str(text).split(';') if x}


def provenance(snapshot, inventory, path):
    key = path.relative_to(snapshot).as_posix()
    row = inventory[key]
    return {'path': str(path), 'remote_source_path': row['source_path'],
            'inventory_bytes': row['bytes'], 'current_bytes': path.stat().st_size,
            'bytes_match': path.stat().st_size == row['bytes']}


def input_audit(root, snapshot):
    constructs = read_constructs(snapshot)
    antibody = read_antibodies(snapshot)
    annotations = read_csv(snapshot / 'data/antibody_imgt_numbering.csv')
    cdr = read_csv(snapshot / 'data/cdr_annotation.csv')
    for chain, label in (('H', 'VH'), ('L', 'VL')):
        rows = sorted((r for r in annotations if r['chain'] == label), key=lambda r: int(r['sequence_index']))
        if ([int(r['sequence_index']) for r in rows] != list(range(1, len(antibody[chain])+1))
                or ''.join(r['amino_acid'] for r in rows) != antibody[chain]):
            raise ValueError(f'{label} snapshot IMGT annotation does not match the frozen antibody')
        covered = []
        for region in (r for r in cdr if r['chain'] == label):
            start, end = int(region['sequence_start']), int(region['sequence_end'])
            if (region['sequence'] != antibody[chain][start-1:end]
                    or int(region['length']) != end-start+1
                    or any(r['region'] != region['region'] for r in rows[start-1:end])):
                raise ValueError(f'{label} snapshot CDR/FR annotation mismatch')
            covered.extend(range(start, end+1))
        if sorted(covered) != list(range(1, len(antibody[chain])+1)):
            raise ValueError(f'{label} snapshot CDR/FR coverage is not complete and unique')
    return {'construct_lengths': {k: len(v) for k, v in constructs.items()},
            'antibody_lengths': {k: len(v) for k, v in antibody.items()}, 'annotations_valid': True,
            'imgt_equals_release': annotations == read_csv(root / 'data/antibody_imgt_numbering.csv'),
            'cdr_equals_release': cdr == read_csv(root / 'data/cdr_annotation.csv')}


def score(value, label):
    if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{label} is missing, nonfinite, or outside 0-1')
    return value


def validate_artifacts(snapshot, inventory, construct, target, directory):
    pdb = directory / f'{target}_model_0.pdb'
    conf_path = directory / f'confidence_{target}_model_0.json'
    pae_path = directory / f'pae_{target}_model_0.npz'
    record_path = directory.parents[1] / 'processed/records' / f'{target}.json'
    sources = {name: provenance(snapshot, inventory, path) for name, path in
               [('pdb', pdb), ('confidence', conf_path), ('pae', pae_path), ('record', record_path)]}
    if not all(row['bytes_match'] for row in sources.values()):
        raise ValueError('A copied artifact size differs from server_source_inventory.json')
    mapping = build_structure_mapping(pdb, construct, snapshot, require_complete=True)
    expected = {'A': CONSTRUCTS[construct][1], 'H': 119, 'L': 106}
    lengths = {chain: mapping[chain]['metadata']['observed_length'] for chain in 'AHL'}
    if lengths != expected:
        raise ValueError(f'Unexpected chain lengths: {lengths}')
    for line in pdb.read_text(encoding='utf-8').splitlines():
        if line.startswith('ATOM  '):
            values = [float(line[x:x+8]) for x in (30, 38, 46)]
            bfactor = float(line[60:66])
            if not all(math.isfinite(v) for v in values) or not math.isfinite(bfactor) or not 0 <= bfactor <= 100:
                raise ValueError('Nonfinite PDB coordinate or pLDDT outside 0-100')
    confidence = json.loads(conf_path.read_text(encoding='utf-8'))
    for field in SCORES:
        score(confidence.get(field), field)
    for field in ('ligand_iptm', 'protein_iptm', 'complex_iplddt'):
        if field in confidence:
            score(confidence[field], field)
    for field in ('complex_pde', 'complex_ipde'):
        value = confidence.get(field)
        if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
            raise ValueError(f'{field} must be finite and nonnegative in angstrom units')
    with np.load(pae_path, allow_pickle=False) as archive:
        pae = archive['pae']
        total = sum(expected.values())
        if (pae.shape != (total, total) or pae.dtype.kind not in 'fiu'
                or not np.isfinite(pae).all() or np.any(pae < 0)):
            raise ValueError(f'Invalid PAE: shape {pae.shape}; expected {(total, total)}')
        pae_range = (float(pae.min()), float(pae.max()))
    record = json.loads(record_path.read_text(encoding='utf-8'))
    chains = record.get('chains', [])
    if record.get('id') != target or len(chains) != 3 or {r['chain_name'] for r in chains} != set('AHL'):
        raise ValueError('Processed record is not this A/H/L target')
    ids = {row['chain_name']: str(row['chain_id']) for row in chains}
    if len(set(ids.values())) != 3:
        raise ValueError('Processed chain IDs are not unique')
    for row in chains:
        if row['num_residues'] != expected[row['chain_name']] or row.get('valid') is not True or row.get('mol_type') != 0:
            raise ValueError('Processed chain identity/length validity differs from the PDB contract')
    for left in 'AHL':
        score(confidence.get('chains_ptm', {}).get(ids[left]), f'{left} chain pTM')
        for right in 'AHL':
            score(confidence.get('pair_chains_iptm', {}).get(ids[left], {}).get(ids[right]), f'{left}-{right} ipTM')
    result = {'sequence_identity': True, 'confidence_valid': True, 'pae_valid': True,
              'chain_A_length': lengths['A'], 'chain_H_length': 119, 'chain_L_length': 106,
              'pae_shape': f'{total}x{total}', 'pae_min_A': pae_range[0], 'pae_max_A': pae_range[1],
              'processed_record_valid': True, 'chain_id_mapping': json.dumps(ids, sort_keys=True),
              'chain_pair_confidence_valid': True, 'source_inventory_bytes_match': True,
              'confidence_path': str(conf_path), 'remote_confidence_path': sources['confidence']['remote_source_path'],
              'pae_path': str(pae_path), 'remote_pae_path': sources['pae']['remote_source_path'],
              'processed_record_path': str(record_path), 'remote_processed_record_path': sources['record']['remote_source_path']}
    return result, confidence, mapping


def pair_key(row):
    return int(row['antigen_local']), row['antibody_chain'], int(row['antibody_local'])


def compare_saved(summary, contacts, cdr, saved, saved_contacts):
    differences = []
    for field in (*SCORES, 'epitope_residues_4p5', 'paratope_residues_4p5',
                  'cdr_paratope_residues_4p5', 'hotspot_hits_4p5', 'hotspot_coverage', 'hotspot_fraction_in_interface'):
        if abs(float(summary[field]) - float(saved[field])) > 1e-12:
            differences.append({'field': field, 'saved': saved[field], 'recomputed': summary[field]})
    if residue_set(summary['epitope_project_residues']) != residue_set(saved['epitope_project_residues']):
        differences.append({'field': 'epitope_project_residues', 'saved': saved['epitope_project_residues'],
                            'recomputed': summary['epitope_project_residues']})
    actual = {pair_key(row): row for row in contacts}
    previous = {pair_key(row): row for row in saved_contacts}
    if len(actual) != len(contacts) or len(previous) != len(saved_contacts):
        raise ValueError('Duplicate residue-pair keys in recalculated or saved contacts')
    absent, added = sorted(previous.keys()-actual.keys()), sorted(actual.keys()-previous.keys())
    field_differences, distance_outliers, rounded_differences, deltas = [], [], [], []
    for key in sorted(actual.keys() & previous.keys()):
        new, old = actual[key], previous[key]
        for field in ('antigen_project', 'antigen_resname', 'antibody_resname', 'antibody_region', 'contact_4p5', 'contact_5p0'):
            if str(new[field]) != str(old[field]):
                field_differences.append({'pair': key, 'field': field, 'saved': old[field], 'recomputed': new[field]})
        delta = abs(new['min_distance_A']-float(old['min_distance_A']))
        deltas.append(delta)
        detail = {'pair': key, 'saved_distance_A': float(old['min_distance_A']),
                  'recomputed_distance_A': new['min_distance_A'], 'absolute_difference_A': delta}
        if delta > DISTANCE_TOLERANCE_A:
            distance_outliers.append(detail)
        if round(new['min_distance_A'], 3) != float(old['min_distance_A']):
            rounded_differences.append(detail)
    old_paratope = {(r['antibody_chain'], int(r['antibody_local'])) for r in saved_contacts if int(r['contact_4p5']) == 1}
    new_paratope = {(r['antibody_chain'], int(r['antibody_local'])) for r in contacts if int(r['contact_4p5']) == 1}
    old_cdr = {(chain, region): len({int(r['antibody_local']) for r in saved_contacts
                                   if int(r['contact_4p5']) == 1 and r['antibody_chain'] == chain and r['antibody_region'] == region})
               for chain in 'HL' for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')}
    cdr_equal = all(row['paratope_residues_4p5'] == old_cdr[(row['antibody_chain'], row['region'])] for row in cdr)
    matching = not any((differences, absent, added, field_differences, distance_outliers)) and new_paratope == old_paratope and cdr_equal
    result = {'saved_summary_matches': not differences, 'summary_difference_count': len(differences),
              'saved_contact_pairs': len(previous), 'recomputed_contact_pairs': len(actual),
              'missing_recomputed_contact_pairs': len(absent), 'added_recomputed_contact_pairs': len(added),
              'contact_field_difference_count': len(field_differences),
              'distance_outside_saved_precision_count': len(distance_outliers),
              'distance_rounded_3dp_difference_count': len(rounded_differences),
              'distance_max_abs_difference_A': max(deltas, default=0.0), 'distance_abs_tolerance_A': DISTANCE_TOLERANCE_A,
              'paratope_set_matches_saved_contacts': new_paratope == old_paratope, 'cdr_counts_match_saved_contacts': cdr_equal,
              'saved_comparison_status': 'matched_at_saved_precision' if matching else 'differences'}
    details = {'summary_differences': differences, 'missing_recomputed_pairs': absent,
               'added_recomputed_pairs': added, 'contact_field_differences': field_differences,
               'distance_outliers': distance_outliers, 'distance_rounding_differences': rounded_differences}
    return result, details


def write_checkpoint(path, records, input_checks, started, complete):
    result = {'started_utc': started, 'updated_utc': now(), 'planned_models': 60,
              'processed_models': len(records), 'audit_complete': complete,
              'raw_artifact_validation_passed': sum(r['raw_validation_status'] == 'passed' for r in records),
              'saved_comparison_matched': sum(r['saved_comparison_status'] == 'matched_at_saved_precision' for r in records),
              'models_with_differences_or_errors': sum(r['status'] != 'passed' for r in records),
              'group_counts': dict(Counter(f"{r['construct']}/{r['condition']}" for r in records)),
              'input_audit': input_checks, 'records': records, 'gpu_tasks_started': 0,
              'verification_method': 'CPU parsing and geometry/contact recalculation; direct comparisons; no digests',
              'distance_comparison_note': 'Old contact distances are stored at 3 decimals; up to 0.00055 A allows that rounding and legacy float32 coordinates. All rounded-value differences are separately retained; threshold flags and residue sets are compared exactly.'}
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--snapshot', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    snapshot = (args.snapshot or root / 'results/full181_supplement/server_baseline_snapshot').resolve()
    out = (args.output or root / 'results/full181_supplement/legacy_raw_audit').resolve()
    if out == snapshot or snapshot in out.parents:
        raise ValueError('Audit outputs must be outside the read-only snapshot')
    out.mkdir(parents=True, exist_ok=True)
    inventory_path = snapshot / 'server_source_inventory.json'
    inventory_rows = json.loads(inventory_path.read_text(encoding='utf-8'))
    inventory = {row['archive_path']: row for row in inventory_rows}
    if len(inventory) != len(inventory_rows):
        raise ValueError('Duplicate archive paths in source inventory')
    inputs = input_audit(root, snapshot)
    prior_path = snapshot / 'results/12_complex_prediction/phase1_prior_deduplicated.csv'
    prior = {int(row['project_residue']) for row in read_csv(prior_path)}
    if not prior:
        raise ValueError('Frozen snapshot Phase1 prior is empty; cannot compare saved prior metrics')
    inputs['frozen_prior'] = provenance(snapshot, inventory, prior_path)
    saved_tables = {}
    for condition, folder in CONDITIONS.items():
        base = snapshot / 'results/13_interface_analysis' / folder
        table = read_csv(base / 'complex_model_summary.csv')
        index = {(r['route'], int(r['seed'])): r for r in table}
        expected_keys = {(route, seed) for route, _ in CONSTRUCTS.values() for seed in range(1, 16)}
        if len(index) != len(table) or set(index) != expected_keys:
            raise ValueError(f'Saved {condition} model summary is not exactly the planned 30 unique models')
        saved_tables[condition] = (base, index, read_csv(base / 'interface_contacts.csv'))
    fields = {
        'legacy_model_summary.csv': list(dict.fromkeys(SUMMARY + EXTRA + ['route', 'hotspot_hits_4p5', 'hotspot_coverage', 'hotspot_fraction_in_interface', 'remote_confidence_path', 'remote_pae_path', 'processed_record_path', 'remote_processed_record_path', 'chain_id_mapping'])),
        'legacy_geometry_qc.csv': GEOMETRY + EXTRA,
        'legacy_cdr_participation.csv': CDR + EXTRA,
        'legacy_interface_contacts.csv': CONTACTS + EXTRA,
        'validation.csv': list(dict.fromkeys(VALIDATION)),
    }
    handles, writers = {}, {}
    started, records = now(), []
    try:
        for name, columns in fields.items():
            handles[name] = (out / name).open('w', encoding='utf-8', newline='')
            writers[name] = csv.DictWriter(handles[name], fieldnames=columns)
            writers[name].writeheader()
        write_checkpoint(out / 'validation.json', records, inputs, started, False)
        for condition, folder in CONDITIONS.items():
            base, summaries, saved_contacts = saved_tables[condition]
            for construct, (route, length) in CONSTRUCTS.items():
                target = route if condition == 'empty' else route + '_msa'
                for seed in range(1, 16):
                    directory = snapshot / 'results/12_complex_prediction' / folder / f'{route}_seed{seed}' / f'boltz_results_{target}' / 'predictions' / target
                    pdb = directory / f'{target}_model_0.pdb'
                    remote = inventory.get(pdb.relative_to(snapshot).as_posix(), {}).get('source_path', 'NA')
                    common = {'construct': construct, 'condition': condition, 'seed': seed, 'source_path': str(pdb),
                              'remote_source_path': remote, 'server_source_path': remote,
                              'source_inventory_path': str(inventory_path), 'source_condition': folder,
                              'raw_validation_status': 'pending', 'saved_comparison_status': 'pending',
                              'evidence_level': 'recomputed_server_original_unrelaxed_PDB'}
                    validation = dict(common, status='pending', failure_reason='')
                    summary, contacts, qc, cdr = None, [], None, []
                    try:
                        checked, confidence, mapping = validate_artifacts(snapshot, inventory, construct, target, directory)
                        validation.update(checked, raw_validation_status='passed')
                        summary, contacts, qc, cdr = analyze_model(snapshot, pdb, condition, seed, Path(checked['confidence_path']), construct)
                        epitope = residue_set(summary['epitope_project_residues'])
                        hits = len(epitope & prior)
                        summary.update(route=route, hotspot_hits_4p5=hits, hotspot_coverage=hits/len(prior),
                                       hotspot_fraction_in_interface=hits/len(epitope) if epitope else 0,
                                       pae_path=checked['pae_path'],
                                       **{key: checked[key] for key in ('remote_confidence_path', 'remote_pae_path', 'processed_record_path', 'remote_processed_record_path', 'chain_id_mapping')})
                        old_contacts = [r for r in saved_contacts if r['route'] == route and int(r['seed']) == seed]
                        comparison, details = compare_saved(summary, contacts, cdr, summaries[(route, seed)], old_contacts)
                        validation.update(comparison, details=details,
                                          status='passed' if comparison['saved_comparison_status'] == 'matched_at_saved_precision' else 'differences_from_saved')
                        for key, path in [('saved_summary', base / 'complex_model_summary.csv'), ('saved_contacts', base / 'interface_contacts.csv')]:
                            origin = provenance(snapshot, inventory, path)
                            validation[key+'_path'] = str(path)
                            validation['remote_'+key+'_path'] = origin['remote_source_path']
                    except Exception as exc:
                        validation.update(status='failed', failure_reason=repr(exc), traceback=traceback.format_exc())
                        if validation['raw_validation_status'] == 'pending':
                            validation['raw_validation_status'] = 'failed'
                        if validation['saved_comparison_status'] == 'pending':
                            validation['saved_comparison_status'] = 'not_completed'
                    extra = {key: validation[key] for key in EXTRA}
                    summary = dict(summary or dict(common, status='failed_validation', failure_reason=validation['failure_reason']), **extra)
                    qc = dict(qc or dict(common, status='failed_validation', failure_reason=validation['failure_reason']), **extra)
                    if not cdr:
                        cdr = [dict(common, status='failed_validation', failure_reason=validation['failure_reason'], antibody_chain=chain, region=region)
                               for chain in 'HL' for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')]
                    products = {'legacy_model_summary.csv': [summary], 'legacy_geometry_qc.csv': [qc],
                                'legacy_cdr_participation.csv': [dict(r, **extra) for r in cdr],
                                'legacy_interface_contacts.csv': [dict(r, **extra) for r in contacts],
                                'validation.csv': [validation]}
                    for name, rows in products.items():
                        writers[name].writerows({key: row.get(key, 'NA') for key in fields[name]} for row in rows)
                        handles[name].flush()
                    records.append(validation)
                    print(f"{len(records)}/60 {construct} {condition} seed{seed}: {validation['status']}", flush=True)
                    if len(records) % 5 == 0:
                        write_checkpoint(out / 'validation.json', records, inputs, started, len(records) == 60)
        write_checkpoint(out / 'validation.json', records, inputs, started, True)
    finally:
        for handle in handles.values():
            handle.close()
    errors = sum(row['status'] != 'passed' for row in records)
    print(json.dumps({'models_audited': len(records), 'passed': len(records)-errors, 'differences_or_errors': errors,
                      'contact_pairs': sum(row.get('recomputed_contact_pairs', 0) for row in records), 'output': str(out)}))
    raise SystemExit(1 if errors else 0)


if __name__ == '__main__':
    main()
