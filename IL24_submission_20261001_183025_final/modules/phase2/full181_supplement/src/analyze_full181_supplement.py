#!/usr/bin/env python3
"""Contact/CDR/confidence/geometry export for accepted Full181 formal models.

Missing tasks get explicit NA records, never artificial empty epitopes.
The optional output directory keeps simultaneous formal seed-1 analyses isolated.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from full181_contract import build_structure_mapping

ROOT = Path(__file__).resolve().parents[1]
NA = 'NA'


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: row.get(k, NA) for k in fields} for row in rows)


def antibody_annotation(root):
    rows = {}
    for row in read_csv(root / 'data/antibody_imgt_numbering.csv'):
        chain = {'VH': 'H', 'VL': 'L'}[row['chain']]
        rows[(chain, int(row['sequence_index']))] = row
    return rows


def protein_atoms(pdb, mapping):
    import numpy as np
    chosen = {}
    for line in pdb.read_text(encoding='utf-8').splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith('ATOM  '):
            continue
        chain, name = line[21], line[12:16].strip()
        element = line[76:78].strip() or name.lstrip('0123456789')[0]
        if chain not in 'AHL' or element.upper() in ('H', 'D'):
            continue
        resid = (' ', int(line[22:26]), line[26] or ' ')
        meta = mapping[chain]['residue_map'][resid]
        occupancy = float(line[54:60].strip() or 0)
        alt = line[16]
        atom = {'chain': chain, 'resid': resid, 'name': name,
                'resname': line[17:20].strip(), 'local': meta['local_position'],
                'project': meta['project_position'], 'amino_acid': meta['amino_acid'],
                'xyz': np.array([float(line[x:x+8]) for x in (30, 38, 46)], dtype=float),
                'bfactor': float(line[60:66]), 'occupancy': occupancy,
                'alt_preference': int(alt in (' ', 'A'))}
        key = (chain, resid, name)
        if key not in chosen or (occupancy, atom['alt_preference']) > (chosen[key]['occupancy'], chosen[key]['alt_preference']):
            chosen[key] = atom
    return list(chosen.values())


def residue_pairs(atoms):
    import numpy as np
    from scipy.spatial import cKDTree
    antigen = [a for a in atoms if a['chain'] == 'A']
    antibody = [a for a in atoms if a['chain'] in 'HL']
    if not antigen or not antibody:
        raise ValueError('No protein heavy atoms on one side of the interface')
    neighbors = cKDTree(np.array([a['xyz'] for a in antibody])).query_ball_point(
        np.array([a['xyz'] for a in antigen]), 5.0)
    pairs = {}
    for i, js in enumerate(neighbors):
        for j in js:
            a, b = antigen[i], antibody[j]
            distance = float(np.linalg.norm(a['xyz'] - b['xyz']))
            if distance >= 5.0:
                continue
            key = (a['local'], b['chain'], b['local'])
            if key not in pairs or distance < pairs[key][2]:
                pairs[key] = (a, b, distance)
    return [pairs[key] for key in sorted(pairs)]


def geometry(atoms):
    import numpy as np
    from scipy.spatial import cKDTree
    by_chain = {c: [a for a in atoms if a['chain'] == c] for c in 'AHL'}
    breaks = cn_outliers = missing = intra = inter = antigen_inter = 0
    for chain, chain_atoms in by_chain.items():
        residues = {}
        for atom in chain_atoms:
            residues.setdefault(atom['local'], {})[atom['name']] = atom
        ordered = [residues[k] for k in sorted(residues)]
        missing += sum(len(set(('N', 'CA', 'C', 'O')) - set(r)) for r in ordered)
        for left, right in zip(ordered, ordered[1:]):
            if 'CA' in left and 'CA' in right:
                breaks += int(np.linalg.norm(left['CA']['xyz'] - right['CA']['xyz']) > 4.5)
            if 'C' in left and 'N' in right:
                distance = np.linalg.norm(left['C']['xyz'] - right['N']['xyz'])
                cn_outliers += int(distance < 1.1 or distance > 1.6)
        tree = cKDTree(np.array([a['xyz'] for a in chain_atoms]))
        intra += sum(abs(chain_atoms[i]['local'] - chain_atoms[j]['local']) > 1
                     for i, j in tree.query_pairs(1.8))
    for i, left in enumerate('AHL'):
        for right in 'AHL'[i+1:]:
            lt = cKDTree(np.array([a['xyz'] for a in by_chain[left]]))
            rt = cKDTree(np.array([a['xyz'] for a in by_chain[right]]))
            count = sum(len(js) for js in lt.query_ball_tree(rt, 1.8))
            inter += count
            if left == 'A':
                antigen_inter += count
    return {'chain_breaks': breaks, 'peptide_cn_outliers': cn_outliers,
            'interchain_clashes_within_1p8': inter,
            'antigen_antibody_clashes_within_1p8': antigen_inter,
            'nonlocal_intrachain_clashes_within_1p8': intra,
            'missing_backbone_atoms': missing,
            'geometry_warning': int(any((breaks, cn_outliers, inter, intra, missing))),
            'thresholds': 'legacy: CA>4.5A; CN<1.1 or >1.6A; cKDTree radius1.8A; nonlocal local-distance>1',
            'qc_source': 'src/assess_relaxed_geometry.py; same thresholds; explicit local mapping',
            'exclusion_policy': 'warnings retained, never excluded after seeing results'}


def chain_pair_confidence(directory, confidence):
    result = {k: NA for k in ('A_H_iptm', 'H_A_iptm', 'A_L_iptm', 'L_A_iptm')}
    result['chain_pair_mapping_status'] = 'unavailable_without_processed_record'
    target = directory.name
    record_path = directory.parents[1] / 'processed/records' / f'{target}.json'
    if not record_path.is_file():
        return result
    record = json.loads(record_path.read_text(encoding='utf-8'))
    chains = record.get('chains', [])
    names = {row.get('chain_name'): str(row.get('chain_id')) for row in chains}
    scores = confidence.get('pair_chains_iptm', confidence.get('chain_pair_iptm', {}))
    if not all(chain in names for chain in 'AHL') or not scores:
        return result
    for left, right in (('A', 'H'), ('H', 'A'), ('A', 'L'), ('L', 'A')):
        value = scores.get(names[left], {}).get(names[right])
        if isinstance(value, (int, float)) and math.isfinite(value):
            result[f'{left}_{right}_iptm'] = value
    result['chain_pair_mapping_status'] = str(record_path)
    return result


SUMMARY = ['construct', 'condition', 'seed', 'source_path', 'status', 'failure_reason',
           'confidence_score', 'ptm', 'iptm', 'complex_plddt', 'confidence_json_scale',
           'pdb_ca_plddt_mean', 'pdb_ca_plddt_min', 'pdb_ca_plddt_max', 'pdb_plddt_scale',
           'epitope_residues_4p5', 'paratope_residues_4p5', 'cdr_paratope_residues_4p5',
           'epitope_project_residues', 'paratope_local_residues', 'epitope_project_residues_5p0',
           'no_contact_model', 'geometry_warning', 'A_H_iptm', 'H_A_iptm', 'A_L_iptm', 'L_A_iptm',
           'chain_pair_mapping_status', 'confidence_path', 'pae_path']
CONTACTS = ['construct', 'condition', 'seed', 'source_path', 'antigen_local', 'antigen_project',
            'antigen_pdb_resseq', 'antigen_pdb_insertion', 'antigen_resname', 'antibody_chain',
            'antibody_local', 'antibody_pdb_resseq', 'antibody_pdb_insertion', 'antibody_resname',
            'antibody_imgt_position', 'antibody_imgt_insertion', 'antibody_region',
            'min_distance_A', 'contact_4p5', 'contact_5p0']
GEOMETRY = ['construct', 'condition', 'seed', 'source_path', 'status', 'failure_reason',
            'chain_breaks', 'peptide_cn_outliers', 'interchain_clashes_within_1p8',
            'antigen_antibody_clashes_within_1p8', 'nonlocal_intrachain_clashes_within_1p8',
            'missing_backbone_atoms', 'geometry_warning', 'thresholds', 'qc_source', 'exclusion_policy']
CDR = ['construct', 'condition', 'seed', 'source_path', 'status', 'failure_reason',
       'antibody_chain', 'region', 'paratope_residues_4p5', 'participates_4p5', 'local_residues']


def analyze_model(root, pdb, condition, seed, confidence_path=None, construct='Full181'):
    mapping = build_structure_mapping(pdb, construct, root, require_complete=True)
    atoms = protein_atoms(pdb, mapping)
    annotation = antibody_annotation(root)
    pairs = residue_pairs(atoms)
    common = dict(construct=construct, condition=condition, seed=seed, source_path=str(pdb))
    contacts, epitope, expanded, paratope = [], set(), set(), set()
    for antigen, antibody, distance in pairs:
        region = annotation[(antibody['chain'], antibody['local'])]
        row = dict(common, antigen_local=antigen['local'], antigen_project=antigen['project'],
                   antigen_pdb_resseq=antigen['resid'][1], antigen_pdb_insertion=antigen['resid'][2].strip(),
                   antigen_resname=antigen['resname'], antibody_chain=antibody['chain'],
                   antibody_local=antibody['local'], antibody_pdb_resseq=antibody['resid'][1],
                   antibody_pdb_insertion=antibody['resid'][2].strip(), antibody_resname=antibody['resname'],
                   antibody_imgt_position=region['imgt_position'], antibody_imgt_insertion=region['imgt_insertion'],
                   antibody_region=region['region'], min_distance_A=distance,
                   contact_4p5=int(distance < 4.5), contact_5p0=int(distance < 5.0))
        contacts.append(row)
        expanded.add(antigen['project'])
        if distance < 4.5:
            epitope.add(antigen['project'])
            paratope.add((antibody['chain'], antibody['local']))
    qc = dict(common, status='complete', failure_reason='', **geometry(atoms))
    ca = [a['bfactor'] for a in atoms if a['name'] == 'CA']
    expected_ca = sum(mapping[c]['metadata']['observed_length'] for c in 'AHL')
    if len(ca) != expected_ca or any(not math.isfinite(v) or not 0 <= v <= 100 for v in ca):
        raise ValueError('CA pLDDT values missing/nonfinite/outside frozen 0-100 PDB scale')
    summary = dict(common, status='complete', failure_reason='',
                   epitope_residues_4p5=len(epitope), paratope_residues_4p5=len(paratope),
                   cdr_paratope_residues_4p5=sum(annotation[k]['region'].startswith('CDR') for k in paratope),
                   epitope_project_residues=';'.join(map(str, sorted(epitope))),
                   epitope_project_residues_5p0=';'.join(map(str, sorted(expanded))),
                   paratope_local_residues=';'.join(f'{c}:{p}' for c, p in sorted(paratope)),
                   no_contact_model=int(not epitope), geometry_warning=qc['geometry_warning'],
                   pdb_ca_plddt_mean=sum(ca)/len(ca), pdb_ca_plddt_min=min(ca),
                   pdb_ca_plddt_max=max(ca), pdb_plddt_scale='0-100', confidence_json_scale='0-1')
    if confidence_path:
        confidence = json.loads(confidence_path.read_text(encoding='utf-8'))
        summary.update({k: confidence.get(k, NA) for k in ('confidence_score', 'ptm', 'iptm', 'complex_plddt')})
        summary.update(chain_pair_confidence(pdb.parent, confidence))
        summary['confidence_path'] = str(confidence_path)
    cdr = []
    for chain in 'HL':
        for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'):
            local = sorted(p for c, p in paratope if c == chain and annotation[(c, p)]['region'] == region)
            cdr.append(dict(common, status='complete', failure_reason='', antibody_chain=chain,
                            region=region, paratope_residues_4p5=len(local), participates_4p5=int(bool(local)),
                            local_residues=';'.join(map(str, local))))
    return summary, contacts, qc, cdr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--condition', choices=['empty', 'formal_msa'])
    parser.add_argument('--seed', type=int, choices=range(1, 16))
    parser.add_argument('--audit-legacy', action='store_true', help='Recalculate only the existing representative PDB subset, read-only')
    args = parser.parse_args()
    out = args.output or ROOT / 'results/full181_supplement'
    out.mkdir(parents=True, exist_ok=True)
    if args.audit_legacy:
        audit = out / 'interface_audit'
        audit.mkdir(parents=True, exist_ok=True)
        inventory = read_csv(ROOT / 'results/full181_supplement/baseline_audit/available_structure_inventory.csv')
        saved = {(r['route'], int(r['seed'])): r for r in read_csv(ROOT / 'results/complex_model_summary.csv')}
        seen, records = set(), []
        for model in inventory:
            key = (model['construct'], model['condition'], int(model['seed']))
            if key in seen:
                continue
            seen.add(key)
            summary, pairs, qc, cdr = analyze_model(ROOT, Path(model['source_path']), key[1], key[2], construct=key[0])
            route = {'Native155': 'native_blind', 'Immunogen134': 'immunogen_blind'}[key[0]]
            old = saved[(route, key[2])]
            agreement = {f'{metric}_matches_saved': str(summary[metric]) == old[metric]
                         for metric in ('epitope_project_residues', 'epitope_residues_4p5',
                                        'paratope_residues_4p5', 'cdr_paratope_residues_4p5')}
            records.append(dict(qc, **agreement, recomputed_contact_pairs_lt5=len(pairs),
                                recalculation_scope='available_original_unrelaxed_PDB_subset_only',
                                actual_geometry_observation=True))
        fields = list(dict.fromkeys(k for r in records for k in r))
        write_csv(audit / 'legacy_raw_subset_interface_qc.csv', records, fields)
        mismatches = sum(not all(v for k, v in r.items() if k.endswith('_matches_saved')) for r in records)
        print(json.dumps({'actual_legacy_PDBs_recalculated': len(records), 'mismatches': mismatches,
                          'old_GPU_tasks_run': 0, 'hashes_computed': 0}))
        if mismatches:
            raise SystemExit(1)
        return
    manifest = read_csv(ROOT / 'results/full181_supplement/task_manifest.csv')
    summaries, contacts, geometries, regions = [], [], [], []
    from run_full181_supplement import validate_artifacts, validate_saved_state, task_dir
    for task in manifest:
        condition, seed = task['condition'], int(task['seed'])
        if args.condition and condition != args.condition or args.seed and seed != args.seed:
            continue
        state_path = task_dir(condition, seed) / 'run_state.json'
        common = dict(construct='Full181', condition=condition, seed=seed, source_path=task['source_path'],
                      status=task['status'], failure_reason=task.get('failure_reason', ''))
        try:
            if not state_path.is_file():
                raise FileNotFoundError('Formal prediction has not started; no Full181 raw artifacts')
            state = json.loads(state_path.read_text(encoding='utf-8'))
            if state.get('status') not in ('complete', 'downstream_blocked'):
                raise ValueError(f'Prediction state is {state.get("status")}')
            artifacts = validate_artifacts(condition, seed)
            validate_saved_state(condition, seed, state, artifacts)
            summary, pairs, qc, cdr = analyze_model(ROOT, Path(artifacts['pdb_path']), condition, seed,
                                                   Path(artifacts['confidence_path']))
            summary['pae_path'] = artifacts['pae_path']
            summaries.append(summary)
            contacts.extend(pairs)
            geometries.append(qc)
            regions.extend(cdr)
        except Exception as exc:
            if task['status'] not in ('blocked_resource', 'planned', 'running'):
                common['status'] = 'failed_analysis_or_validation'
            common['failure_reason'] = task.get('failure_reason') or repr(exc)
            summaries.append(dict(common))
            geometries.append(dict(common))
            for chain in 'HL':
                for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'):
                    regions.append(dict(common, antibody_chain=chain, region=region))
    write_csv(out / 'full181_model_summary.csv', summaries, SUMMARY)
    write_csv(out / 'full181_interface_contacts.csv', contacts, CONTACTS)
    write_csv(out / 'full181_geometry_qc.csv', geometries, GEOMETRY)
    write_csv(out / 'full181_cdr_participation.csv', regions, CDR)
    print(json.dumps({'planned_rows': len(summaries), 'completed_observations': sum(r['status'] == 'complete' for r in summaries),
                      'contact_rows': len(contacts), 'output': str(out)}, ensure_ascii=True))
    if any(row['status'] == 'failed_analysis_or_validation' for row in summaries):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
