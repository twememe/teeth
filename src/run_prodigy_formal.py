#!/usr/bin/env python3
"""Run PRODIGY on formal antigen(A)-antibody(H,L) complexes."""
import argparse, csv, os, re, shutil, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUTPUT_ROOT=Path(os.environ.get("PHASE2_OUTPUT_ROOT",ROOT/"results")).resolve()
ap=argparse.ArgumentParser(); ap.add_argument("--set",choices=["formal","formal_msa"],default="formal"); args=ap.parse_args()
base=OUTPUT_ROOT/"12_complex_prediction"/args.set
out=OUTPUT_ROOT/"15_affinity"; out.mkdir(parents=True,exist_ok=True)
prodigy=shutil.which("prodigy")
if not prodigy: raise SystemExit("prodigy-prot 2.4.0 console command not found")
rows=[]
for pdb in sorted(base.glob("**/*_model_0.pdb")):
    run=pdb.parts[pdb.parts.index(args.set)+1]; route,seed=run.rsplit("_seed",1)
    p=subprocess.run([prodigy,str(pdb),"--selection","A","H,L"],text=True,capture_output=True,check=True)
    dg=re.search(r"binding affinity .*?:\s+(-?[0-9.]+)",p.stdout)
    kd=re.search(r"dissociation constant .*?:\s+([0-9.eE+-]+)",p.stdout)
    nc=re.search(r"No\. of intermolecular contacts[:]\s+(\d+)",p.stdout)
    rows.append({"route":route,"seed":seed,"prodigy_dG_kcal_mol":dg.group(1),"prodigy_Kd_M_25C":kd.group(1),"intermolecular_contacts":nc.group(1),"method_note":"relative ranking only; not experimental KD"})
with open(out/f"prodigy_{args.set}.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(f"PRODIGY completed for {len(rows)} models")
