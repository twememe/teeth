#!/usr/bin/env python3
"""Strict, dependency-light audit for the compact release."""
import csv, hashlib, re, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
FIELDS={'candidate_id','track','task','route','seed','msa_condition','sequence_or_structure','source_pdb','structure_status','relaxation_status','clash_qc','model_version','ipTM','complex_pLDDT','PRODIGY_dG_kcal_mol','hotspot_coverage','ranking_basis','notes'}
REQUIRED=['README.md','requirements.txt','requirements-full.txt','environment.yml','MODEL_CARD.md','REPRODUCIBILITY.md','THIRD_PARTY_NOTICES.md','predict.py','run.sh','run_full.sh','results.csv','data/README.md','data/msa/msa_manifest.json','data/phase1_prior_deduplicated.csv','models/MODEL_ASSETS.csv','scripts/fetch_model_assets.py','logs/LOG_INDEX.csv','logs/README.md','docs/DATA_PROVENANCE.csv','docs/MODEL_PROVENANCE.csv','docs/SOFTWARE_VERSIONS.csv','docs/STRUCTURE_FILES.md','docs/README_reproduce.md','licenses/Boltz-LICENSE.txt','licenses/ImmuneBuilder-LICENSE.txt','src/build_recomputed_results.py']
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def fail(message): raise SystemExit(message)
def main():
    missing=[x for x in REQUIRED if not (ROOT/x).is_file()]
    if missing: fail('Missing required files: '+', '.join(missing))
    first=subprocess.run(['bash','run.sh'],cwd=ROOT,capture_output=True,text=True)
    if first.returncode: fail('first run.sh failed: '+first.stderr)
    h1=sha(ROOT/'results.csv')
    second=subprocess.run(['bash','run.sh'],cwd=ROOT,capture_output=True,text=True)
    if second.returncode: fail('second run.sh failed: '+second.stderr)
    if h1!=sha(ROOT/'results.csv'): fail('run.sh is not deterministic')
    rows=list(csv.DictReader((ROOT/'results.csv').open(encoding='utf-8')))
    if len(rows)!=10 or not rows: fail(f'Expected 10 candidates, found {len(rows)}')
    if not FIELDS.issubset(rows[0]): fail('Missing result fields: '+', '.join(sorted(FIELDS-set(rows[0]))))
    if any(r['track']!='赛道一：AI大分子与多肽药物设计' or not r['task'] or not r['ranking_basis'] for r in rows): fail('Invalid track, task, or ranking basis')
    if len({r['candidate_id'] for r in rows})!=10 or len({r['source_pdb'] for r in rows})!=10: fail('IDs or structures are not unique')
    for r in rows:
        p=ROOT/r['source_pdb']
        if not p.is_file() or p.resolve().is_relative_to(ROOT.resolve()) is False: fail(f'Invalid structure path: {r["source_pdb"]}')
        chains={line[21].strip() for line in p.read_text(errors='ignore').splitlines() if line.startswith('ATOM')}
        if not {'A','H','L'}.issubset(chains): fail(f'Unexpected chains in {p}: {chains}')
        for k in ('seed','ipTM','complex_pLDDT','PRODIGY_dG_kcal_mol','hotspot_coverage'): float(r[k])
    if (ROOT/'results_immunogen.csv').exists() or (ROOT/'results_native.csv').exists(): fail('Obsolete route-split result file present')
    req=(ROOT/'requirements.txt').read_text(); full=(ROOT/'requirements-full.txt').read_text()
    if 'prodigy-prot==2.4.0' not in req or 'PRODIGY==' in req+full or 'OpenMM==8.6.0' not in req: fail('Dependency records are inconsistent')
    assets=list(csv.DictReader((ROOT/'models/MODEL_ASSETS.csv').open(encoding='utf-8')))
    if len(assets)!=6 or any(not r['sha256'] or len(r['sha256'])!=64 or not r['primary_url'] or int(r['size_bytes'])<=0 for r in assets): fail('Invalid MODEL_ASSETS.csv')
    for r in assets[:2]:
        p=ROOT/'models'/r['group']/r['asset']
        if p.exists() and (p.stat().st_size!=int(r['size_bytes']) or sha(p)!=r['sha256']): fail(f'Invalid bundled model asset: {p}')
    logs=list(csv.DictReader((ROOT/'logs/LOG_INDEX.csv').open(encoding='utf-8')))
    if len(logs)!=73: fail(f'Expected 73 indexed logs, found {len(logs)}')
    if any(not (ROOT/'logs'/r['log_file']).is_file() for r in logs): fail('Log index references a missing file')
    names={r['log_file'] for r in logs}
    expected={f'boltz_formal_msa_{route}_seed{i}.log' for route in ('native_blind','immunogen_blind') for i in range(1,16)}
    expected|={f'boltz_formal_{route}_seed{i}.log' for route in ('native_blind','immunogen_blind') for i in range(1,16)}
    expected|={f'boltz_glyco_{route}_seed{i}.log' for route in ('native_blind','immunogen_blind') for i in range(1,4)}
    expected|={f'boltz_msa_{route}_seed{i}.log' for route in ('native_blind','immunogen_blind') for i in range(1,4)}|{'formal_msa_driver.log'}
    if names!=expected: fail('Log coverage or naming is incomplete')
    provenance=list(csv.DictReader((ROOT/'docs/DATA_PROVENANCE.csv').open(encoding='utf-8')))
    required_cols={'file','source','version','acquired_date','purpose','used_for_training','used_as_prediction_input','posterior_only','license','sha256','notes'}
    if not provenance or not required_cols.issubset(provenance[0]): fail('DATA_PROVENANCE.csv fields incomplete')
    data_files={str(p.relative_to(ROOT)).replace('\\','/') for p in (ROOT/'data').rglob('*') if p.is_file()}
    provenance_files={r['file'] for r in provenance}
    if data_files!=provenance_files: fail('Data provenance coverage mismatch: '+str(sorted(data_files^provenance_files)))
    for r in provenance:
        p=ROOT/r['file']
        if not all(r[c] for c in required_cols) or sha(p)!=r['sha256']: fail('Incomplete or stale provenance row: '+r['file'])
    forbidden=re.compile(r'/(?:home|mnt)/|[A-Za-z]:\\\\')
    scan=[ROOT/'README.md',ROOT/'docs',ROOT/'models',ROOT/'scripts',ROOT/'src',ROOT/'run.sh',ROOT/'run_full.sh']
    for base in scan:
        files=[base] if base.is_file() else [p for p in base.rglob('*') if p.is_file()]
        for p in files:
            if forbidden.search(p.read_text(errors='ignore')): fail('Absolute path found: '+str(p.relative_to(ROOT)))
    for p in ROOT.rglob('*'):
        if p.name=='__pycache__' or p.suffix=='.pyc' or p.name.endswith(('.partial','.lock')): fail('Temporary artifact found: '+str(p.relative_to(ROOT)))
    isolation=['run_imgt_numbering.py','step11_antibody_structure.py','generate_boltz_inputs.py','run_boltz_formal_msa.sh','analyze_boltz_complexes.py','run_prodigy_formal.py','rank_formal_msa_models.py','summarize_formal_msa.py','build_recomputed_results.py']
    for name in isolation:
        if 'PHASE2_OUTPUT_ROOT' not in (ROOT/'src'/name).read_text(errors='ignore'): fail('Output isolation missing in src/'+name)
    test_out=ROOT/'work/test_preflight_output'
    check=subprocess.run(['bash','run_full.sh','--preflight-only','--output-dir','work/test_preflight_output'],cwd=ROOT,capture_output=True,text=True)
    if check.returncode not in (0,1) or 'OUTPUT_TARGET=work/test_preflight_output' not in check.stdout or test_out.exists(): fail('Preflight output isolation test failed')
    failure=subprocess.run(['bash','run_full.sh','--self-test-failure'],cwd=ROOT,capture_output=True,text=True)
    if failure.returncode==0: fail('Full entry failed to propagate a child failure')
    manifest_lines=(ROOT/'RELEASE_SHA256SUMS.txt').read_text().splitlines(); listed=set()
    for line in manifest_lines:
        h,rel=line.split('  ',1); rel=rel.removeprefix('./'); listed.add(rel); p=ROOT/rel
        if not p.is_file() or sha(p)!=h: fail('Checksum failure: '+rel)
    actual={str(p.relative_to(ROOT)).replace('\\','/') for p in ROOT.rglob('*') if p.is_file() and p.name!='RELEASE_SHA256SUMS.txt'}
    if listed!=actual: fail('SHA-256 manifest coverage mismatch: '+str(sorted(listed^actual)))
    print(f'VALIDATION_OK candidates=10 logs=73 assets=6 data_files={len(data_files)} deterministic_sha256={h1}')
if __name__=='__main__': main()
