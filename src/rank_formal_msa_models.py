#!/usr/bin/env python3
"""Rank formal-MSA models and select confidence/consensus representatives."""

from __future__ import annotations

import csv, os
import shutil
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
OUTPUT_ROOT=Path(os.environ.get("PHASE2_OUTPUT_ROOT",ROOT/"results")).resolve()
summary=list(csv.DictReader((OUTPUT_ROOT/"13_interface_analysis/formal_msa/complex_model_summary.csv").open()))
aff={(r["route"],r["seed"]):r for r in csv.DictReader((OUTPUT_ROOT/"15_affinity/prodigy_formal_msa.csv").open())}
out=OUTPUT_ROOT/"14_quality"
reps=out/"representatives_msa";reps.mkdir(parents=True,exist_ok=True)


def epitope(row):
    return {int(x) for x in row["epitope_project_residues"].split(";") if x}


rows=[]
for route in sorted({r["route"] for r in summary}):
    selected=[r for r in summary if r["route"]==route]
    for row in selected:
        a=epitope(row)
        scores=[]
        for other in selected:
            if other is row: continue
            b=epitope(other);scores.append(len(a&b)/len(a|b) if a|b else 1.0)
        row["mean_epitope_jaccard_to_others"]=sum(scores)/len(scores)
    confidence=sorted(selected,key=lambda r:float(r["confidence_score"]),reverse=True)
    consensus=sorted(selected,key=lambda r:(float(r["mean_epitope_jaccard_to_others"]),float(r["confidence_score"])),reverse=True)
    informed=sorted(selected,key=lambda r:(float(r["hotspot_coverage"]),float(r["confidence_score"])),reverse=True)
    for row in selected:
        iptm=float(row["iptm"]); quality="high" if iptm>=.6 else "medium" if iptm>=.45 else "low"
        rows.append({"route":route,"seed":row["seed"],"quality_class":quality,
                     "confidence_score":row["confidence_score"],"iptm":row["iptm"],
                     "complex_plddt":row["complex_plddt"],"hotspot_coverage":row["hotspot_coverage"],
                     "mean_epitope_jaccard_to_others":row["mean_epitope_jaccard_to_others"],
                     "prodigy_dG_kcal_mol":aff[(route,row["seed"])]["prodigy_dG_kcal_mol"],
                     "confidence_rank":confidence.index(row)+1,"consensus_rank":consensus.index(row)+1,
                     "phase1_informed_rank":informed.index(row)+1})
    choices={"confidence":confidence[0],"consensus":consensus[0]}
    if route=="native_blind": choices["phase1_informed"]=informed[0]
    for label,row in choices.items():
        source=next((OUTPUT_ROOT/"12_complex_prediction/formal_msa"/f"{route}_seed{row['seed']}").glob("**/*_model_0.pdb"))
        shutil.copy2(source,reps/f"{route}_msa_{label}_representative_seed{row['seed']}.pdb")

with (out/"model_ranking_formal_msa.csv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(f"Ranked {len(rows)} models; wrote {len(list(reps.glob('*.pdb')))} representatives")
