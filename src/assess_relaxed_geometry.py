#!/usr/bin/env python3
"""Verify geometry and coordinate drift after restrained minimization."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree


ROOT=Path(__file__).resolve().parents[1]
OUTPUT_ROOT=Path(os.environ.get("PHASE2_OUTPUT_ROOT",ROOT/"results")).resolve()
RELAXED=OUTPUT_ROOT/"14_quality/relaxed"
parser=PDBParser(QUIET=True)


def get_model(path): return next(parser.get_structure(path.stem,path).get_models())


def heavy_atoms(model,chain): return [a for a in model[chain].get_atoms() if a.element!="H"]


def geometry(model):
    breaks=cn_outliers=0
    for chain in ("A","H","L"):
        residues=[r for r in model[chain] if r.id[0]==" "]
        for left,right in zip(residues,residues[1:]):
            if "CA" in left and "CA" in right and np.linalg.norm(left["CA"].coord-right["CA"].coord)>4.5: breaks+=1
            if "C" in left and "N" in right:
                distance=np.linalg.norm(left["C"].coord-right["N"].coord)
                if distance<1.1 or distance>1.6: cn_outliers+=1
    inter=0
    for index,first in enumerate(("A","H","L")):
        for second in ("A","H","L")[index+1:]:
            a=np.asarray([x.coord for x in heavy_atoms(model,first)])
            b=np.asarray([x.coord for x in heavy_atoms(model,second)])
            inter+=sum(len(x) for x in cKDTree(a).query_ball_tree(cKDTree(b),1.8))
    intra=0
    for chain in ("A","H","L"):
        atoms=heavy_atoms(model,chain);coords=np.asarray([a.coord for a in atoms])
        for left,right in cKDTree(coords).query_pairs(1.8):
            if abs(atoms[left].get_parent().id[1]-atoms[right].get_parent().id[1])>1: intra+=1
    return breaks,cn_outliers,inter,intra


def backbone(model):
    return np.asarray([a.coord for c in ("A","H","L") for r in model[c] if r.id[0]==" " for name in ("N","CA","C","O") if name in r for a in (r[name],)])


def main():
    rows=[]
    for report in sorted(RELAXED.glob("*_relaxation.json")):
        metadata=json.loads(report.read_text())
        source=Path(metadata["source"])
        destination=Path(metadata["destination"])
        before=get_model(source);after=get_model(destination)
        before_bb=backbone(before);after_bb=backbone(after)
        count=min(len(before_bb),len(after_bb))
        drift=float(np.sqrt(np.mean(np.sum((after_bb[:count]-before_bb[:count])**2,axis=1))))
        b0,c0,i0,n0=geometry(before);b1,c1,i1,n1=geometry(after)
        rows.append({"model":destination.stem,"backbone_atom_count":count,"backbone_coordinate_rmsd_A":drift,
                     "before_chain_breaks":b0,"after_chain_breaks":b1,"before_peptide_cn_outliers":c0,"after_peptide_cn_outliers":c1,
                     "before_interchain_clashes_lt1p8":i0,"after_interchain_clashes_lt1p8":i1,
                     "before_nonlocal_intrachain_clashes_lt1p8":n0,"after_nonlocal_intrachain_clashes_lt1p8":n1,
                     "initial_potential_energy_kj_mol":metadata["initial_potential_energy_kj_mol"],
                     "final_potential_energy_kj_mol":metadata["final_potential_energy_kj_mol"]})
    with (OUTPUT_ROOT/"14_quality/relaxed_geometry_qc.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(f"Assessed {len(rows)} relaxed representatives")


if __name__=="__main__": main()
