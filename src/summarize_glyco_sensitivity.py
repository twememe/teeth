#!/usr/bin/env python3
"""Compare matched formal-MSA models with and without one N74-linked NAG."""

from __future__ import annotations

import csv, os
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree


ROOT=Path(__file__).resolve().parents[1]
OUTPUT_ROOT=Path(os.environ.get("PHASE2_OUTPUT_ROOT",ROOT/"results")).resolve()
INTERFACE=OUTPUT_ROOT/"13_interface_analysis"
COMPLEX=OUTPUT_ROOT/"12_complex_prediction"
OUT=OUTPUT_ROOT/"14_quality"
parser=PDBParser(QUIET=True)


def read(condition):
    return list(csv.DictReader((INTERFACE/condition/"complex_model_summary.csv").open()))


def epitope(row):
    return {int(x) for x in row["epitope_project_residues"].split(";") if x}


def nag_distance(route,seed):
    pdb=next((COMPLEX/"glyco_sensitivity"/f"{route}_seed{seed}").glob("**/*_model_0.pdb"))
    model=next(parser.get_structure("x",pdb).get_models())
    nag=np.asarray([a.coord for a in model["G"].get_atoms() if a.element!="H"])
    ab=np.asarray([a.coord for c in ("H","L") for a in model[c].get_atoms() if a.element!="H"])
    return float(min(cKDTree(ab).query(x,k=1)[0] for x in nag))


def main():
    base=read("formal_msa")
    glyco=read("glyco_sensitivity")
    rows=[]
    for route in ("native_blind","immunogen_blind"):
        for seed in ("1","2","3"):
            unmodified=next(r for r in base if r["route"]==route and r["seed"]==seed)
            modified=next(r for r in glyco if r["route"]==route and r["seed"]==seed)
            a,b=epitope(unmodified),epitope(modified)
            rows.append({
                "route":route,"seed":seed,
                "unmodified_iptm":unmodified["iptm"],"n74_nag_iptm":modified["iptm"],
                "delta_iptm":float(modified["iptm"])-float(unmodified["iptm"]),
                "unmodified_complex_plddt":unmodified["complex_plddt"],"n74_nag_complex_plddt":modified["complex_plddt"],
                "delta_complex_plddt":float(modified["complex_plddt"])-float(unmodified["complex_plddt"]),
                "unmodified_hotspot_coverage":unmodified["hotspot_coverage"],"n74_nag_hotspot_coverage":modified["hotspot_coverage"],
                "delta_hotspot_coverage":float(modified["hotspot_coverage"])-float(unmodified["hotspot_coverage"]),
                "paired_epitope_jaccard":len(a&b)/len(a|b) if a|b else 1.0,
                "unmodified_n74_contact":int(74 in a),"n74_nag_n74_contact":int(74 in b),
                "nag_min_heavy_atom_distance_to_antibody_A":nag_distance(route,seed),
            })
    with (OUT/"glyco_sensitivity_paired.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=[]
    for route in ("native_blind","immunogen_blind"):
        selected=[r for r in rows if r["route"]==route]
        result={"route":route,"n_pairs":len(selected)}
        for metric in ("delta_iptm","delta_complex_plddt","delta_hotspot_coverage","paired_epitope_jaccard","nag_min_heavy_atom_distance_to_antibody_A"):
            values=np.asarray([float(r[metric]) for r in selected])
            result[f"{metric}_mean"]=float(values.mean());result[f"{metric}_sd"]=float(values.std(ddof=1))
        result["unmodified_n74_contact_models"]=sum(r["unmodified_n74_contact"] for r in selected)
        result["n74_nag_n74_contact_models"]=sum(r["n74_nag_n74_contact"] for r in selected)
        summary.append(result)
    with (OUT/"glyco_sensitivity_summary.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    print("Wrote paired N74-NAG sensitivity analysis")


if __name__=="__main__": main()
