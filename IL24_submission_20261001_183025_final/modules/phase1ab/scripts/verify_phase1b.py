#!/usr/bin/env python3
"""Fresh final verification and freeze manifest for IL-24 Phase 1B."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scripts.candidates_phase1b import audit_phase1a_manifest
    from scripts.integration import read_fasta
except ModuleNotFoundError:
    from candidates_phase1b import audit_phase1a_manifest
    from integration import read_fasta


ROOT = Path(__file__).resolve().parents[1]
STOP_FILES = [
    "results/05b_discotope_restored/discotope_residue_scores.csv",
    "results/07b_integration_with_3d/residue_feature_table_3d.csv",
    "results/08b_candidates_with_3d/candidate_regions_phase1b_preunblind.csv",
    "results/08b_candidates_with_3d/candidate_stability_matrix_phase1b.csv",
    "results/08b_candidates_with_3d/phase1a_vs_phase1b_candidates.csv",
    "reports/phase1b_preunblind_report.md",
]
FIGURES = [f"results/08b_candidates_with_3d/0{i}_{name}.png" for i, name in [
    (1, "all_evidence_tracks_phase1b"),
    (2, "candidate_profile_phase1b"),
    (3, "robustness_heatmap_phase1b"),
    (4, "bepipred_vs_discotope"),
    (5, "phase1a_vs_phase1b_profile"),
    (6, "phase1a_vs_phase1b_candidate_rank"),
]]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str, checks: list[str]) -> None:
    if not condition:
        raise AssertionError(message)
    checks.append(message)


def finite(frame: pd.DataFrame, columns: list[str]) -> bool:
    return np.isfinite(frame[columns].astype(float).to_numpy()).all()


def pdb_bfactor_map(path: Path) -> tuple[dict[int, float], str]:
    values: dict[int, float] = {}
    remark = ""
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("REMARK 900"):
            remark += line + "\n"
        if line.startswith("ATOM  ") and line[21:22] == "A":
            position = int(line[22:26])
            values.setdefault(position, float(line[60:66]))
    return values, remark


def freeze_paths() -> list[Path]:
    paths = [
        ROOT / "protocol_phase1b_amendment.yaml",
        ROOT / "reports/phase1b_amendment.md",
        ROOT / "reports/step05b_discotope_restored.md",
        ROOT / "reports/step07b_integration_with_3d.md",
        ROOT / "reports/phase1b_preunblind_report.md",
    ]
    for directory in (
        ROOT / "results/phase1a_baseline",
        ROOT / "results/05b_discotope_restored",
        ROOT / "results/07b_integration_with_3d",
        ROOT / "results/08b_candidates_with_3d",
    ):
        paths.extend(path for path in directory.rglob("*") if path.is_file())
    for pattern in (
        "install_dt3_phase1b_attempt*.log",
        "phase1b_network_preflight.json",
        "phase1b_model_load_smoke.log",
        "phase1b_esm_if1_range_download.log",
        "discotope_phase1b_gpu*.log",
        "discotope_phase1b_gpu_usage.csv",
    ):
        paths.extend(path for path in (ROOT / "logs").glob(pattern) if path.is_file())
    return sorted(set(paths), key=lambda path: path.relative_to(ROOT).as_posix())


def main() -> None:
    checks: list[str] = []
    reference = read_fasta(ROOT / "inputs/reference_181.fasta")
    require(len(reference) == 181, "reference length is exactly 181", checks)
    for relative in STOP_FILES:
        require((ROOT / relative).is_file(), f"stop file exists: {relative}", checks)

    dt = pd.read_csv(ROOT / STOP_FILES[0])
    require(len(dt) == 181, "Step 5B table has 181 rows", checks)
    require(dt["position"].tolist() == list(range(1, 182)), "Step 5B numbering is 1-181", checks)
    require("".join(dt["AA"].astype(str)) == reference, "Step 5B AA sequence is exact", checks)
    require(finite(dt, ["DiscoTope_raw", "DiscoTope_calibrated", "DiscoTope_RSA", "DiscoTope_pLDDT"]), "Step 5B continuous fields are finite", checks)

    table = pd.read_csv(ROOT / STOP_FILES[1])
    require(len(table) == 181, "Step 7B table has 181 rows", checks)
    require(table["position"].tolist() == list(range(1, 182)), "Step 7B numbering is 1-181", checks)
    require("".join(table["AA"].astype(str)) == reference, "Step 7B AA sequence is exact", checks)
    numeric = [
        "pLDDT", "BepiPred_raw", "BepiPred_percentile", "DiscoTope_raw",
        "DiscoTope_percentile", "EpitopeConsensus", "SASA", "RSA_raw",
        "RSA_percentile", "Conservation_raw", "Conservation_percentile",
        "FinalScore_Phase1A", "FinalScore_Phase1B", "FinalScore_4way_sensitivity",
    ]
    require(finite(table, numeric), "Step 7B numeric fields are finite", checks)
    e = (table["BepiPred_percentile"] + table["DiscoTope_percentile"]) / 2
    final_b = (e + table["RSA_percentile"] + table["Conservation_percentile"]) / 3
    four_way = (table["BepiPred_percentile"] + table["DiscoTope_percentile"] + table["RSA_percentile"] + table["Conservation_percentile"]) / 4
    require(np.allclose(table["EpitopeConsensus"], e, atol=1e-12), "EpitopeConsensus formula is exact", checks)
    require(np.allclose(table["FinalScore_Phase1B"], final_b, atol=1e-12), "Phase 1B primary formula is exact", checks)
    require(np.allclose(table["FinalScore_4way_sensitivity"], four_way, atol=1e-12), "four-way sensitivity formula is exact", checks)
    phase1a_table = pd.read_csv(ROOT / "results/phase1a_baseline/residue_feature_table.csv")
    error = float((table["FinalScore_Phase1A"] - phase1a_table["FinalScore"]).abs().max())
    require(error < 1e-12, "BepiPred-only reproduces frozen Phase 1A", checks)

    expected_window_rows = {15: 167, 20: 162, 25: 157}
    for scale, count in expected_window_rows.items():
        frame = pd.read_csv(ROOT / f"results/08b_candidates_with_3d/window_{scale}_scores_3d.csv")
        require(len(frame) == count, f"{scale}-aa table retains all {count} windows", checks)
        require((frame["length"] == scale).all(), f"{scale}-aa table uses fixed window length", checks)
        require(finite(frame, ["mean_score", "FinalScore_Phase1B_mean"]), f"{scale}-aa window scores are finite", checks)

    candidates = pd.read_csv(ROOT / STOP_FILES[2])
    require(len(candidates) == 5, "Phase 1B candidate table has Top 5 families", checks)
    candidate_required = [
        "candidate_id", "start", "end", "length", "FinalScore_Phase1B",
        "BepiPred_mean", "DiscoTope_mean", "EpitopeConsensus_mean", "RSA_mean",
        "Conservation_mean", "median_pLDDT", "fraction_pLDDT_ge_70",
        "scale_stability", "LOO_stability",
    ]
    require(set(candidate_required).issubset(candidates.columns), "candidate schema contains all required fields", checks)
    require((candidates["bootstrap_repetitions"] == 100).all(), "candidate conservation bootstrap is exactly 100", checks)
    stability = pd.read_csv(ROOT / STOP_FILES[3])
    require(len(stability) == 20, "LOO matrix has 5 candidates x 4 settings", checks)
    require(set(stability["setting"]) == {"Full", "-Epitope", "-Surface", "-Conservation"}, "LOO settings are grouped evidence settings", checks)
    sensitivity = pd.read_csv(ROOT / "results/08b_candidates_with_3d/epitope_model_sensitivity.csv")
    require(len(sensitivity) == 15, "epitope model sensitivity retains 3 x Top 5 families", checks)
    require(set(sensitivity["model_setting"]) == {"Full", "BepiPred-only", "DiscoTope-only"}, "epitope model sensitivity settings are complete", checks)

    for relative in FIGURES:
        path = ROOT / relative
        require(path.is_file() and path.stat().st_size > 50_000, f"figure is nonempty: {relative}", checks)

    final_pdb = ROOT / "results/08b_candidates_with_3d/IL24_181_phase1b_candidate_score.pdb"
    dt_pdb = ROOT / "results/08b_candidates_with_3d/IL24_181_discotope_score.pdb"
    final_values, final_remark = pdb_bfactor_map(final_pdb)
    dt_values, dt_remark = pdb_bfactor_map(dt_pdb)
    require(set(final_values) == set(range(1, 182)), "Phase 1B FinalScore PDB maps all 181 residues", checks)
    require(set(dt_values) == set(range(1, 182)), "DiscoTope PDB maps all 181 residues", checks)
    require("PHASE1B FINALSCORE" in final_remark, "Phase 1B PDB labels B-factor semantics", checks)
    require("DISCOTOPE PERCENTILE" in dt_remark, "DiscoTope PDB labels B-factor semantics", checks)
    require(max(abs(final_values[pos] / 100 - table.loc[pos - 1, "FinalScore_Phase1B"]) for pos in final_values) < 0.0001, "Phase 1B PDB B-factors match scores", checks)
    require(max(abs(dt_values[pos] / 100 - table.loc[pos - 1, "DiscoTope_percentile"]) for pos in dt_values) < 0.0001, "DiscoTope PDB B-factors match percentiles", checks)

    report = (ROOT / STOP_FILES[-1]).read_text(encoding="utf-8")
    require("Step 9 was not performed" in report, "final report explicitly prohibits Step 9", checks)
    require("Historical experimental intervals remain blinded" in report, "final report confirms historical blinding", checks)

    phase1a_audit = audit_phase1a_manifest(ROOT, ROOT / "logs/phase1_preunblind_freeze_sha256.csv")
    require(phase1a_audit["all_unchanged"] and phase1a_audit["entries"] == 68, "Phase 1A freeze replay is 68/68 unchanged", checks)

    manifest_path = ROOT / "logs/phase1b_preunblind_freeze_sha256.csv"
    frozen_paths = freeze_paths()
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        for path in frozen_paths:
            writer.writerow({
                "path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })

    payload = {
        "status": "PASS",
        "checks_passed": len(checks),
        "checks": checks,
        "phase1a_audit": phase1a_audit,
        "phase1b_freeze_manifest": str(manifest_path.relative_to(ROOT)),
        "phase1b_freeze_entries": len(frozen_paths),
        "phase1b_freeze_manifest_sha256": sha256(manifest_path),
        "stop_files": STOP_FILES,
        "step9_performed": False,
    }
    (ROOT / "logs/phase1b_final_verification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
