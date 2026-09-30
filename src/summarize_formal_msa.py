#!/usr/bin/env python3
"""Summarize and compare the 15-seed formal empty-MSA and MSA ensembles."""

from __future__ import annotations

import csv, os
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = Path(os.environ.get("PHASE2_OUTPUT_ROOT", ROOT / "results")).resolve()
INTERFACE = OUTPUT_ROOT / "13_interface_analysis"
OUT = OUTPUT_ROOT / "14_quality"


def read(condition: str) -> list[dict[str, str]]:
    return list(csv.DictReader((INTERFACE / condition / "complex_model_summary.csv").open(encoding="utf-8")))


def epitope(row: dict[str, str]) -> set[int]:
    return {int(x) for x in row["epitope_project_residues"].split(";") if x}


def jaccards(rows: list[dict[str, str]]) -> list[float]:
    values=[]
    for left,right in combinations(rows,2):
        a,b=epitope(left),epitope(right)
        values.append(len(a&b)/len(a|b) if a|b else 1.0)
    return values


def main() -> None:
    rows=[]
    conditions=[x for x in ("formal","formal_msa") if (INTERFACE/x/"complex_model_summary.csv").is_file()]
    if not conditions: raise FileNotFoundError("no formal interface summaries found")
    for condition in conditions:
        data=read(condition)
        condition_consensus={}
        for route in ("native_blind","immunogen_blind"):
            selected=[row for row in data if row["route"]==route]
            counts=Counter(x for row in selected for x in epitope(row))
            consensus={x for x,n in counts.items() if n>=len(selected)/2}
            condition_consensus[route]=consensus
            jac=np.asarray(jaccards(selected))
            result={"condition":condition,"route":route,"n_models":len(selected),
                    "pairwise_epitope_jaccard_mean":float(jac.mean()),
                    "pairwise_epitope_jaccard_sd":float(jac.std(ddof=1)),
                    "epitope_consensus_ge50pct":";".join(map(str,sorted(consensus)))}
            for metric in ("iptm","complex_plddt","hotspot_coverage","epitope_residues_4p5"):
                values=np.asarray([float(row[metric]) for row in selected])
                result[f"{metric}_mean"]=float(values.mean())
                result[f"{metric}_sd"]=float(values.std(ddof=1))
            rows.append(result)
        shared=set(range(27,161))
        a=condition_consensus["native_blind"]&shared
        b=condition_consensus["immunogen_blind"]&shared
        cross=len(a&b)/len(a|b) if a|b else 1.0
        for row in rows:
            if row["condition"]==condition:
                row["native_immunogen_consensus_jaccard"]=cross

    fields=[]
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/"formal_msa_comparison.csv").open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    print("Wrote formal 15-seed MSA comparison")


if __name__ == "__main__":
    main()
