#!/usr/bin/env python3
"""Close Phase 2 structural-QC gaps using the formal Boltz ensembles."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
COMPLEX = ROOT / "results" / "12_complex_prediction"
QUALITY = ROOT / "results" / "14_quality"
INTERFACE = ROOT / "results" / "13_interface_analysis"
PARSER = PDBParser(QUIET=True)


def models_in(condition: str, route: str) -> list[tuple[int, Path]]:
    result = []
    for path in (COMPLEX / condition).glob(f"{route}_seed*/**/*_model_0.pdb"):
        match = re.search(r"_seed(\d+)", path.as_posix())
        if match:
            result.append((int(match.group(1)), path))
    return sorted(result)


def model(path: Path):
    return next(PARSER.get_structure(path.stem, path).get_models())


def atoms_for(structure, chain: str, residues: range, names: tuple[str, ...]) -> np.ndarray:
    coords = []
    for resid in residues:
        residue = structure[chain][(" ", resid, " ")]
        for name in names:
            coords.append(residue[name].coord)
    return np.asarray(coords, dtype=float)


def superpose(reference: np.ndarray, mobile: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    sup = SVDSuperimposer()
    sup.set(reference, mobile)
    sup.run()
    rotation, translation = sup.get_rotran()
    return rotation, translation, float(sup.get_rms())


def transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation + translation


def antibody_ca(structure) -> np.ndarray:
    coords = []
    for chain in ("H", "L"):
        for residue in structure[chain].get_residues():
            if residue.id[0] == " " and "CA" in residue:
                coords.append(residue["CA"].coord)
    return np.asarray(coords, dtype=float)


def pose_comparison() -> None:
    rows = []
    for condition in ("formal", "formal_msa", "msa_sensitivity"):
        natives = models_in(condition, "native_blind")
        immunogens = models_in(condition, "immunogen_blind")
        for native_seed, native_path in natives:
            native = model(native_path)
            native_ag = atoms_for(native, "A", range(1, 135), ("N", "CA", "C"))
            native_ab = antibody_ca(native)
            for immunogen_seed, immunogen_path in immunogens:
                immunogen = model(immunogen_path)
                immunogen_ag = atoms_for(immunogen, "A", range(1, 135), ("N", "CA", "C"))
                immunogen_ab = antibody_ca(immunogen)
                rotation, translation, antigen_rmsd = superpose(native_ag, immunogen_ag)
                moved_ab = transform(immunogen_ab, rotation, translation)
                pose_rmsd = float(np.sqrt(np.mean(np.sum((moved_ab - native_ab) ** 2, axis=1))))
                centroid_shift = float(np.linalg.norm(moved_ab.mean(0) - native_ab.mean(0)))
                _, _, internal_rmsd = superpose(native_ab, immunogen_ab)
                rows.append(
                    {
                        "condition": condition,
                        "native_seed": native_seed,
                        "immunogen_seed": immunogen_seed,
                        "shared_antigen_backbone_atoms": len(native_ag),
                        "shared_antigen_backbone_rmsd_A": round(antigen_rmsd, 4),
                        "antibody_ca_atoms": len(native_ab),
                        "antibody_pose_ca_rmsd_A": round(pose_rmsd, 4),
                        "antibody_centroid_shift_A": round(centroid_shift, 4),
                        "antibody_internal_ca_rmsd_A": round(internal_rmsd, 4),
                        "is_representative_pair": int(condition == "formal" and native_seed == 14 and immunogen_seed == 1),
                    }
                )
    QUALITY.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (QUALITY / "native_immunogen_pose_comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    metrics = [
        "shared_antigen_backbone_rmsd_A",
        "antibody_pose_ca_rmsd_A",
        "antibody_centroid_shift_A",
        "antibody_internal_ca_rmsd_A",
    ]
    summary = []
    for condition in ("formal", "formal_msa", "msa_sensitivity"):
        selected = [row for row in rows if row["condition"] == condition]
        for metric in metrics:
            values = np.asarray([row[metric] for row in selected], dtype=float)
            representatives = [row[metric] for row in selected if row["is_representative_pair"]]
            summary.append(
                {
                    "condition": condition,
                    "metric": metric,
                    "n_pairs": len(values),
                    "mean": round(float(values.mean()), 4),
                    "sd": round(float(values.std(ddof=1)), 4),
                    "median": round(float(np.median(values)), 4),
                    "min": round(float(values.min()), 4),
                    "max": round(float(values.max()), 4),
                    "representative_pair": representatives[0] if representatives else "",
                }
            )
    with (QUALITY / "native_immunogen_pose_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)


def glycosylation_analysis() -> None:
    seq = "".join(
        line.strip()
        for line in (ROOT / "data" / "IL24_full_1_181.fasta").read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    motifs = [(idx + 1, seq[idx : idx + 3]) for idx in range(len(seq) - 2) if seq[idx] == "N" and seq[idx + 1] != "P" and seq[idx + 2] in "ST"]
    rows = []
    for condition in ("formal", "formal_msa", "msa_sensitivity"):
      for route in ("native_blind", "immunogen_blind"):
        for seed, path in models_in(condition, route):
            structure = model(path)
            antibody = np.asarray(
                [atom.coord for chain in ("H", "L") for atom in structure[chain].get_atoms() if atom.element != "H"],
                dtype=float,
            )
            tree = cKDTree(antibody)
            for project_residue, motif in motifs:
                local = project_residue - 26
                if local < 1 or local not in structure["A"]:
                    continue
                residue_coords = np.asarray([atom.coord for atom in structure["A"][(" ", local, " ")] if atom.element != "H"])
                # LOGIC-CLOSEOUT 20260930: this is the distance from the PROTEIN N74 residue to the
                # antibody.  It is NOT a glycan distance, and the field name
                # `glycan_reach_sensitivity_lt15` below is a retained misnomer.  The true
                # NAG-to-antibody distance is computed in summarize_glyco_sensitivity.py.
                distance = float(min(tree.query(coord, k=1)[0] for coord in residue_coords))
                rows.append(
                    {
                        "condition": condition,
                        "route": route,
                        "seed": seed,
                        "project_residue": project_residue,
                        "model_local_residue": local,
                        "sequon": motif,
                        "min_heavy_atom_distance_to_antibody_A": round(distance, 3),
                        "direct_contact_lt4p5": int(distance < 4.5),
                        "near_interface_lt10": int(distance < 10.0),
                        "glycan_reach_sensitivity_lt15": int(distance < 15.0),
                    }
                )
    with (INTERFACE / "glycosylation_interface_analysis.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def geometry_qc() -> None:
    rows = []
    for condition in ("formal", "formal_msa", "msa_sensitivity"):
      for route in ("native_blind", "immunogen_blind"):
        for seed, path in models_in(condition, route):
            structure = model(path)
            chain_breaks = 0
            peptide_bond_outliers = 0
            for chain in ("A", "H", "L"):
                residues = [residue for residue in structure[chain].get_residues() if residue.id[0] == " "]
                for left, right in zip(residues, residues[1:]):
                    if "CA" in left and "CA" in right and np.linalg.norm(left["CA"].coord - right["CA"].coord) > 4.5:
                        chain_breaks += 1
                    if "C" in left and "N" in right:
                        cn = np.linalg.norm(left["C"].coord - right["N"].coord)
                        if cn < 1.1 or cn > 1.6:
                            peptide_bond_outliers += 1

            interchain_clashes = 0
            chain_atoms = {}
            for chain in ("A", "H", "L"):
                chain_atoms[chain] = np.asarray([atom.coord for atom in structure[chain].get_atoms() if atom.element != "H"])
            for idx, first in enumerate(("A", "H", "L")):
                for second in ("A", "H", "L")[idx + 1 :]:
                    interchain_clashes += sum(len(hits) for hits in cKDTree(chain_atoms[first]).query_ball_tree(cKDTree(chain_atoms[second]), 1.8))

            nonlocal_intrachain_clashes = 0
            for chain in ("A", "H", "L"):
                atoms = [atom for atom in structure[chain].get_atoms() if atom.element != "H"]
                coords = np.asarray([atom.coord for atom in atoms])
                for left, right in cKDTree(coords).query_pairs(1.8):
                    res_left = atoms[left].get_parent().id[1]
                    res_right = atoms[right].get_parent().id[1]
                    if abs(res_left - res_right) > 1:
                        nonlocal_intrachain_clashes += 1
            rows.append(
                {
                    "condition": condition,
                    "route": route,
                    "seed": seed,
                    "chain_breaks_ca_gt4p5": chain_breaks,
                    "peptide_cn_outliers_outside_1p1_1p6": peptide_bond_outliers,
                    "interchain_heavy_atom_clashes_lt1p8": interchain_clashes,
                    "nonlocal_intrachain_heavy_atom_clashes_lt1p8": nonlocal_intrachain_clashes,
                    "geometry_pass": int(chain_breaks == 0 and peptide_bond_outliers == 0 and interchain_clashes == 0 and nonlocal_intrachain_clashes == 0),
                }
            )
    with (QUALITY / "structural_geometry_qc.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    pose_comparison()
    glycosylation_analysis()
    geometry_qc()
    print("Wrote pose, glycosylation, and geometry QC outputs.")


if __name__ == "__main__":
    main()
