#!/usr/bin/env python3
"""Analyze Boltz IL-24/antibody complexes using 4.5 and 5.0 A contacts."""
import argparse, csv, json, os
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = Path(os.environ.get("PHASE2_OUTPUT_ROOT", ROOT / "results")).resolve()
ap=argparse.ArgumentParser(); ap.add_argument("--set",choices=["pilot","formal","msa_sensitivity","formal_msa","glyco_sensitivity"],default="pilot"); args=ap.parse_args()
PILOT = OUTPUT_ROOT / "12_complex_prediction" / args.set
OUT = OUTPUT_ROOT / "13_interface_analysis" / args.set
OUT.mkdir(parents=True, exist_ok=True)

prior = set()
with open(ROOT / "data" / "phase1_prior_deduplicated.csv") as f:
    for row in csv.DictReader(f): prior.add(int(row["project_residue"]))

cdr = {}
numbering_cdr = OUTPUT_ROOT / "10_numbering" / "cdr_annotation.csv"
with open(numbering_cdr if numbering_cdr.is_file() else ROOT / "data" / "cdr_annotation.csv") as f:
    for row in csv.DictReader(f):
        chain = "H" if row["chain"] == "VH" else "L"
        for i in range(int(row["sequence_start"]), int(row["sequence_end"])+1):
            cdr[(chain, i)] = row["region"]

summaries=[]; contacts=[]
parser=PDBParser(QUIET=True)
for pdb in sorted(PILOT.glob("**/*_model_0.pdb")):
    route_seed=pdb.parts[pdb.parts.index(args.set)+1]
    route,seed=route_seed.rsplit("_seed",1)
    structure=parser.get_structure("x",pdb)
    model=next(structure.get_models())
    ag_atoms=[]; ag_meta=[]; ab_atoms=[]; ab_meta=[]
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ": continue
            for atom in residue:
                if atom.element == "H": continue
                meta=(chain.id,residue.id[1],residue.resname)
                if chain.id == "A":
                    ag_atoms.append(atom.coord); ag_meta.append(meta)
                elif chain.id in ("H", "L"):
                    ab_atoms.append(atom.coord); ab_meta.append(meta)
    tree=cKDTree(np.asarray(ab_atoms)); pairs=tree.query_ball_point(np.asarray(ag_atoms),5.0)
    pair_min={}
    for i,js in enumerate(pairs):
        for j in js:
            key=(ag_meta[i],ab_meta[j]); d=float(np.linalg.norm(np.asarray(ag_atoms[i])-np.asarray(ab_atoms[j])))
            pair_min[key]=min(pair_min.get(key,999),d)
    for (a,b),d in pair_min.items():
        contacts.append({"route":route,"seed":seed,"antigen_local":a[1],"antigen_project":a[1]+26,"antigen_resname":a[2],"antibody_chain":b[0],"antibody_local":b[1],"antibody_resname":b[2],"antibody_region":cdr.get((b[0],b[1]),"NA"),"min_distance_A":round(d,3),"contact_4p5":int(d<4.5),"contact_5p0":1})
    epi45={a[0][1]+26 for a,d in pair_min.items() if d<4.5}; para45={(a[1][0],a[1][1]) for a,d in pair_min.items() if d<4.5}
    conf=json.loads(next(pdb.parent.glob("confidence_*.json")).read_text())
    hit=epi45 & prior
    summaries.append({"route":route,"seed":seed,"confidence_score":conf["confidence_score"],"ptm":conf["ptm"],"iptm":conf["iptm"],"complex_plddt":conf["complex_plddt"],"epitope_residues_4p5":len(epi45),"paratope_residues_4p5":len(para45),"cdr_paratope_residues_4p5":sum(cdr.get(x,"").startswith("CDR") for x in para45),"hotspot_hits_4p5":len(hit),"hotspot_coverage":len(hit)/len(prior),"hotspot_fraction_in_interface":len(hit)/len(epi45) if epi45 else 0,"epitope_project_residues":";".join(map(str,sorted(epi45)))})

for name,rows in [("interface_contacts.csv",contacts),("complex_model_summary.csv",summaries)]:
    with open(OUT/name,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"Analyzed {len(summaries)} complexes and {len(contacts)} residue pairs")
