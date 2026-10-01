import csv, json, math, re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
EXT=ROOT/'data/external_materials'
OUT=ROOT/'results/16_mechanism'
OUT.mkdir(exist_ok=True)

def atoms(path, chains=None):
    out=[]
    for line in path.read_text(errors='ignore').splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        ch=line[21].strip()
        if chains and ch not in chains: continue
        try: out.append((ch,int(line[22:26]),line[12:16].strip(),tuple(float(line[i:i+8]) for i in (30,38,46))))
        except: pass
    return out

def contacts(pdb, ligand, partners, cutoff=5.0):
    a=atoms(pdb,{ligand}); b=atoms(pdb,set(partners)); c={x:set() for x in partners}
    for ch,r,an,xyz in a:
        for ch2,r2,an2,xyz2 in b:
            if sum((xyz[i]-xyz2[i])**2 for i in range(3)) <= cutoff*cutoff: c.setdefault(ch2,set()).add(r)
    return {k:sorted(v) for k,v in c.items()}

def model_epitope(pdb):
    a=atoms(pdb,{'A'}); h=atoms(pdb,{'H','L'}); s=set()
    for ch,r,an,xyz in a:
        if any(sum((xyz[i]-q[i])**2 for i in range(3))<=25 for _,_,_,q in h): s.add(r)
    return sorted(s)

six=EXT/'6DF3.pdb'
rc=contacts(six,'C',['L','H'])
# PDB chain C uses local residue numbers 1-155; project numbering starts at 27.
rc_project={k:[r+26 for r in v] for k,v in rc.items()}
(OUT/'6DF3_receptor_contacts.json').write_text(json.dumps({'pdb':'6DF3','ligand_chain':'C','receptor_chains':{'L':'IL-22R1','H':'IL-20R2'},'cutoff_A':5.0,'numbering':'project = PDB local + 26','contacts_pdb_local':rc,'contacts_project':rc_project},indent=2),encoding='utf-8')

# 6DF3 chain C is UniProt 52-206; project mature numbering is 27-181.
t198={'literature_residue':'T198','structure_chain_C_auth_residue':198,'project_numbering':173,'mapping_basis':'6DF3/UniProt mature numbering 52-206 to project mature numbering 27-181 (offset -25)','status':'sequence-offset mapping; verify against exact project sequence'}
(OUT/'T198_project_mapping.json').write_text(json.dumps(t198,indent=2),encoding='utf-8')

rows=[]
for p in sorted((ROOT/'results/14_quality/representatives_msa').glob('*.pdb')):
    ep=set(model_epitope(p)); rec=set(rc_project['L'])|set(rc_project['H']); inter=sorted(ep&rec); union=ep|rec
    rows.append({'model':p.name,'epitope_n':len(ep),'receptor_contact_n':len(rec),'overlap_n':len(inter),'jaccard':round(len(inter)/len(union),4) if union else 0,'overlap_residues':';'.join(map(str,inter))})
with (OUT/'6DF3_model_overlap.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)

elisa=list(csv.DictReader((EXT/'ELISA_competition_data.csv').open(encoding='utf-8-sig')))
for r in elisa:
    r['antibody_minus_IgG_pct']=round(float(r['1A6-13-8 抑制率 (%)'])-float(r['同亚型 IgG 抑制率 (%)']),3)
(OUT/'ELISA_external_consistency.json').write_text(json.dumps({'source':'data/external_materials/ELISA_competition_data.csv','rows':elisa,'interpretation':'concentration-dependent inhibition with low isotype control; assay metadata incomplete'},ensure_ascii=False,indent=2),encoding='utf-8')

md=f'''# Phase 2 external-structure analysis\n\n- 6DF3 receptor contacts were computed from the downloaded coordinates using a 5.0 Å heavy-atom distance cutoff. IL-24 is chain C; IL-22R1 is chain L; IL-20R2 is chain H.\n- Literature T198 maps provisionally to project residue T173 by the mature-chain offset (52–206 -> 27–181). Exact sequence confirmation remains required.\n- Model overlap results are in `6DF3_model_overlap.csv`; this is a structural hypothesis, not proof of receptor competition.\n- ELISA data were preserved unchanged and summarized with antibody-minus-isotype inhibition in `ELISA_external_consistency.json`.\n- 6DF3 does not establish the IL-20R1/IL-20R2 complex.\n'''
(OUT/'external_structure_analysis.md').write_text(md,encoding='utf-8')
print('Wrote 6DF3 contacts, model overlap, T198 mapping, ELISA summary and report')
