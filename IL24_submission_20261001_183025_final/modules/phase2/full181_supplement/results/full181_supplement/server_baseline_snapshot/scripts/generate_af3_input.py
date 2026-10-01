#!/usr/bin/env python3
"""Generate auditable Step 12 inputs from frozen FASTA files."""
import hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; DATA = ROOT / "data"; OUT = ROOT / "results" / "12_complex_prediction"; OUT.mkdir(parents=True, exist_ok=True)
def fasta(name):
    seq = ""; header = ""
    for line in (DATA/name).read_text().splitlines():
        if line.startswith(">") : header=line[1:].strip()
        else: seq += line.strip()
    return header, seq
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
_, native = fasta("IL24_native_27_181.fasta"); _, immunogen = fasta("IL24_immunogen_27_160.fasta")
records={}; key=None
for line in (DATA/"VH_VL_clean.fasta").read_text().splitlines():
    if line.startswith(">"):
        key="VH" if "VH" in line else "VL"; records[key]=""
    elif key: records[key]+=line.strip()
hotspots=[("C1-B",52,66),("C2-B",126,140),("C3-B",60,74),("C4-B",143,157),("C5-B",117,131)]
prior={}
for region,start,end in hotspots:
    for r in range(start,end+1): prior.setdefault(r,[]).append(region)
def make(name, antigen, guided=False):
    obj={"name":name,"modelSeeds":list(range(1,21)),"sequences":[{"protein":{"name":"IL24","sequence":antigen}},{"protein":{"name":"VH","sequence":records["VH"]}},{"protein":{"name":"VL","sequence":records["VL"]}}],"metadata":{"numbering":"project 1-181; local mature index = project-26","guided_mode":"Phase1-informed post-ranking" if guided else "blind"}}
    if guided: obj["phase1_prior"]=[{"project_residue":r,"local_residue":r-26,"regions":prior[r],"weight":1.0} for r in sorted(prior)]
    return obj
routes={"12A_native_blind":make("IL24_native_27_181_blind",native),"12B_native_phase1_informed":make("IL24_native_27_181_phase1_informed",native,True),"12C_immunogen_blind":make("IL24_immunogen_27_160_blind",immunogen)}
for route,obj in routes.items():
    d=OUT/route; d.mkdir(parents=True,exist_ok=True); (d/(route+".json")).write_text(json.dumps(obj,indent=2)+"\n")
(OUT/"step12_input_manifest_sha256.csv").write_text("file,sha256\n"+"\n".join(f'{n},{sha(DATA/n)}' for n in ["IL24_native_27_181.fasta","IL24_immunogen_27_160.fasta","VH_VL_clean.fasta"])+"\n")
(OUT/"phase1_prior_deduplicated.csv").write_text("project_residue,local_mature_index,regions,weight\n"+"\n".join(f'{r},{r-26},{";".join(prior[r])},1.0' for r in sorted(prior))+"\n")
print(f"Generated {len(routes)} routes; native={len(native)}, immunogen={len(immunogen)}, VH={len(records['VH'])}, VL={len(records['VL'])}, unique_prior={len(prior)}")
