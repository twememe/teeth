#!/usr/bin/env python3
"""Rank formal models without collapsing distinct evidence into one opaque score."""
import csv, shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
summary=list(csv.DictReader(open(ROOT/"results/13_interface_analysis/formal/complex_model_summary.csv")))
aff={(r["route"],r["seed"]):r for r in csv.DictReader(open(ROOT/"results/15_affinity/prodigy_formal.csv"))}
out=ROOT/"results/14_quality"; reps=out/"representatives"; reps.mkdir(parents=True,exist_ok=True)
rows=[]
for route in sorted({r["route"] for r in summary}):
    rr=[r for r in summary if r["route"]==route]
    blind=sorted(rr,key=lambda r:float(r["confidence_score"]),reverse=True)
    informed=sorted(rr,key=lambda r:(float(r["hotspot_coverage"]),float(r["confidence_score"])),reverse=True)
    energy=sorted(rr,key=lambda r:float(aff[(route,r["seed"])]["prodigy_dG_kcal_mol"]))
    for r in rr:
        iptm=float(r["iptm"]); q="high" if iptm>=.6 else "medium" if iptm>=.45 else "low"
        rows.append({"route":route,"seed":r["seed"],"quality_class":q,"confidence_score":r["confidence_score"],"iptm":r["iptm"],"complex_plddt":r["complex_plddt"],"hotspot_coverage":r["hotspot_coverage"],"prodigy_dG_kcal_mol":aff[(route,r["seed"])]["prodigy_dG_kcal_mol"],"blind_confidence_rank":blind.index(r)+1,"phase1_informed_rank":informed.index(r)+1,"prodigy_rank":energy.index(r)+1})
    choices={"blind_representative":blind[0]}
    if route=="native_blind": choices["phase1_informed_representative"]=informed[0]
    for label,r in choices.items():
        src=next((ROOT/"results/12_complex_prediction/formal"/f"{route}_seed{r['seed']}").glob("**/*_model_0.pdb"))
        shutil.copy2(src,reps/f"{route}_{label}_seed{r['seed']}.pdb")
with open(out/"model_ranking.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(f"Ranked {len(rows)} models; representatives={len(list(reps.glob('*.pdb')))}")
