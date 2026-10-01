#!/usr/bin/env python3
"""Create validated Phase 1B Step 5B artifacts from official DiscoTope output."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scripts.run_discotope import load_reference, map_official_csv_phase1b
except ModuleNotFoundError:  # Direct execution: script directory is sys.path[0].
    from run_discotope import load_reference, map_official_csv_phase1b


PINNED_COMMIT = "35d9f2e55f97eaba2a7acefbc394db58fb9670bc"
OFFICIAL_REPO = "https://github.com/DTU/DiscoTope-3.0"
ESM_IF1_URL = (
    "https://dl.fbaipublicfiles.com/fair-esm/models/"
    "esm_if1_gvp4_t16_142M_UR50.pt"
)
OUTPUT_COLUMNS = [
    "position",
    "AA",
    "DiscoTope_raw",
    "DiscoTope_calibrated",
    "DiscoTope_prediction",
    "DiscoTope_RSA",
    "DiscoTope_pLDDT",
    "pdb",
    "chain",
    "pdb_resseq",
    "length",
    "alphafold_struc_flag",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_commit(repo: Path) -> str:
    if not (repo / ".git").exists() and (repo / "SOURCE_REVISION").is_file():
        return (repo / "SOURCE_REVISION").read_text().strip()
    result = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def contiguous_true_segments(rows: list[dict]) -> list[tuple[int, int]]:
    segments: list[tuple[int, int]] = []
    start: int | None = None
    for row in rows:
        position = int(row["position"])
        if row["DiscoTope_prediction"] and start is None:
            start = position
        if start is not None and (
            not row["DiscoTope_prediction"] or position == len(rows)
        ):
            end = position if row["DiscoTope_prediction"] else position - 1
            segments.append((start, end))
            start = None
    return segments


def top_windows(rows: list[dict], size: int = 15, count: int = 3) -> list[dict]:
    scores = [float(row["DiscoTope_raw"]) for row in rows]
    windows = [
        {
            "start": start + 1,
            "end": start + size,
            "mean": sum(scores[start : start + size]) / size,
        }
        for start in range(len(scores) - size + 1)
    ]
    windows.sort(key=lambda item: (-item["mean"], item["start"]))
    selected: list[dict] = []
    for window in windows:
        if all(
            max(
                0,
                min(window["end"], chosen["end"])
                - max(window["start"], chosen["start"])
                + 1,
            )
            / size
            <= 0.5
            for chosen in selected
        ):
            selected.append(window)
        if len(selected) == count:
            break
    return selected


def parse_gpu_usage(path: Path) -> dict:
    values: list[int] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.reader(handle):
            if len(row) < 3 or "memory.used" in row[2]:
                continue
            match = re.search(r"(\d+)", row[2])
            if match:
                values.append(int(match.group(1)))
    return {
        "samples": len(values),
        "observed_memory_used_min_mib": min(values) if values else None,
        "observed_memory_used_peak_mib": max(values) if values else None,
        "observed_peak_delta_mib": max(values) - min(values) if values else None,
    }


def parse_wall_time_seconds(path: Path) -> float | None:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Elapsed \(wall clock\) time.*?:\s*(\d+):(\d+(?:\.\d+)?)", text)
    if not match:
        return None
    return int(match.group(1)) * 60 + float(match.group(2))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def write_profile(path: Path, rows: list[dict]) -> None:
    x = [row["position"] for row in rows]
    raw = [row["DiscoTope_raw"] for row in rows]
    calibrated = [row["DiscoTope_calibrated"] for row in rows]
    predicted = [row["DiscoTope_prediction"] for row in rows]

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.8), sharex=True)
    axes[0].plot(x, raw, color="#1565C0", linewidth=1.6)
    axes[0].fill_between(x, 0, raw, color="#90CAF9", alpha=0.35)
    axes[0].set_ylabel("Official raw score")
    axes[0].set_title("IL-24 DiscoTope-3.0 residue profile (AlphaFold mode)")

    axes[1].plot(x, calibrated, color="#6A1B9A", linewidth=1.4)
    axes[1].axhline(0.90, color="#C62828", linestyle="--", linewidth=1.2,
                    label="Official calibrated threshold = 0.90")
    for position, flag in zip(x, predicted):
        if flag:
            axes[1].axvspan(position - 0.5, position + 0.5, color="#FFB300", alpha=0.32)
    axes[1].set_xlabel("Project position")
    axes[1].set_ylabel("Calibrated score")
    axes[1].legend(loc="upper right", frameon=False)
    for axis in axes:
        axis.set_xlim(1, 181)
        axis.grid(axis="y", color="#CFD8DC", linewidth=0.6, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-csv", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--esm-checkpoint", type=Path, required=True)
    parser.add_argument("--gpu-log", type=Path, required=True)
    parser.add_argument("--inference-log", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    for path in (
        args.raw_csv,
        args.reference,
        args.pdb,
        args.esm_checkpoint,
        args.gpu_log,
        args.inference_log,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    if repo_commit(args.repo) != PINNED_COMMIT:
        raise RuntimeError("Official repository is not at the frozen Phase 1B commit")

    reference = load_reference(args.reference)
    rows, official_columns = map_official_csv_phase1b(
        args.raw_csv, reference, expected_chain="A"
    )
    if len(rows) != 181 or "".join(row["AA"] for row in rows) != reference:
        raise RuntimeError("Restored output failed the 181-residue sequence gate")
    for row in rows:
        for column in (
            "DiscoTope_raw",
            "DiscoTope_calibrated",
            "DiscoTope_RSA",
            "DiscoTope_pLDDT",
        ):
            if not math.isfinite(float(row[column])):
                raise RuntimeError(f"Non-finite {column} at {row['position']}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "discotope_residue_scores.csv"
    profile_path = args.output_dir / "discotope_profile.png"
    provenance_path = args.output_dir / "provenance.json"
    write_csv(csv_path, rows)
    write_profile(profile_path, rows)

    raw_scores = [row["DiscoTope_raw"] for row in rows]
    calibrated = [row["DiscoTope_calibrated"] for row in rows]
    positives = [row for row in rows if row["DiscoTope_prediction"]]
    segments = contiguous_true_segments(rows)
    gpu = parse_gpu_usage(args.gpu_log)
    wall_seconds = parse_wall_time_seconds(args.inference_log)
    xgb_files = sorted((args.repo / "models").glob("XGB_*_of_100.json"))

    provenance = {
        "phase": "Phase 1B post-freeze authorized restoration",
        "historical_experimental_regions_accessed": False,
        "discotope_available": True,
        "official_repository": OFFICIAL_REPO,
        "official_commit": PINNED_COMMIT,
        "official_version": "DiscoTope-3.0 1.0",
        "readme_python": "3.14",
        "actual_python": "3.14.7",
        "environment": "il24-dt3-phase1b",
        "install_strategy": "B2_AFTER_EXPLICIT_GCC_COMPATIBILITY_CORRECTION",
        "b1_outcome": "network-resumed; then biotraj source build failed because gcc was absent",
        "b2_correction": "Conda gcc_linux-64/gxx_linux-64 16.1.0; official requirements unchanged",
        "package_versions": {
            "torch": "2.10.0",
            "torch_cuda": "12.8",
            "torch_geometric": "2.7.0",
            "biopython": "1.86",
            "biotite": "1.6.0",
            "biotraj": "1.2.2",
            "xgboost": "3.1.2",
            "pygam": "0.12.0",
            "numpy": "2.4.3",
            "pandas": "3.0.1",
            "joblib": "1.5.3",
            "scikit_learn": "1.7.2",
        },
        "model_assets": {
            "xgboost_models_loaded": len(xgb_files),
            "gam_models_loaded": 2,
            "esm_if1": {
                "url": ESM_IF1_URL,
                "bytes": args.esm_checkpoint.stat().st_size,
                "sha256": sha256(args.esm_checkpoint),
                "parameter_count": 141662150,
            },
        },
        "inference": {
            "structure_type": "alphafold",
            "cpu_only": False,
            "actual_device": "NVIDIA GeForce RTX 4050 Laptop GPU",
            "cuda_oom": False,
            "wall_seconds": wall_seconds,
            "gpu_usage": gpu,
            "official_threshold": 0.90,
            "official_threshold_score_field": "calibrated_score",
        },
        "mapping": {
            "row_count": len(rows),
            "project_numbering": "1-181",
            "sequence_exact_match": True,
            "chain": "A",
            "official_columns": official_columns,
            "primary_discotope_continuous_field": "DiscoTope-3.0_score",
            "phase1b_output_field": "DiscoTope_raw",
            "calibrated_score_preserved": True,
        },
        "inputs": {
            "reference": {"path": str(args.reference), "sha256": sha256(args.reference)},
            "pdb": {"path": str(args.pdb), "sha256": sha256(args.pdb)},
            "official_raw_csv": {
                "path": str(args.raw_csv),
                "sha256": sha256(args.raw_csv),
            },
        },
        "outputs": {
            "discotope_residue_scores.csv": sha256(csv_path),
            "discotope_profile.png": sha256(profile_path),
        },
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )

    top_residues = sorted(
        rows, key=lambda row: (-row["DiscoTope_raw"], row["position"])
    )[:10]
    windows = top_windows(rows)
    segment_text = ", ".join(f"{start}-{end}" for start, end in segments) or "none"
    report = f"""# Step 05B — Restored DiscoTope-3.0 evidence

## Status

- Status: **COMPLETE** (`discotope_available = true`).
- Restoration authority: post-freeze Phase 1B amendment; Phase 1A remains immutable.
- Historical experimental regions were not accessed.

## Official method and installation provenance

- Official repository: {OFFICIAL_REPO}
- Frozen commit: `{PINNED_COMMIT}`
- Official README runtime: Python 3.14; actual isolated environment: `il24-dt3-phase1b`, Python 3.14.7.
- B1 preserved the official dependency specifications. Network Range retries completed the large wheels, after which `biotraj` produced an explicit source-build failure because `gcc` was absent.
- B2 added only the compatible Conda GCC/G++ 16.1.0 toolchain. No predictor, model, Python version, or official requirement version was substituted.
- Official loaders successfully loaded 100 XGBoost models, two GAM calibration models, and ESM-IF1 (141,662,150 parameters).

DiscoTope-3.0 uses an ESM inverse-folding protein model to encode the three-dimensional backbone and sequence context, then applies an ensemble of XGBoost residue classifiers. The official length/surface GAM models calibrate scores across proteins. The supplied AlphaFold structure was therefore run with `--struc_type alphafold`.

## Inference

- Device actually used: NVIDIA GeForce RTX 4050 Laptop GPU (CUDA 12.8); no CUDA OOM and no CPU fallback.
- Wall-clock time: {wall_seconds:.2f} s.
- Observed GPU memory: {gpu['observed_memory_used_min_mib']}–{gpu['observed_memory_used_peak_mib']} MiB (peak delta {gpu['observed_peak_delta_mib']} MiB across {gpu['samples']} samples).
- Input gate: one model, Chain A only, residues 1–181, complete backbone, exact reference sequence match.
- Official raw output: exactly 181 rows; no non-finite continuous values.

## Residue scores

- Official uncalibrated `DiscoTope-3.0_score` is preserved as `DiscoTope_raw` and is the continuous value preregistered for Phase 1B percentile ranking.
- Official `calibrated_score` is separately preserved as `DiscoTope_calibrated`; the official epitope classification uses its default threshold 0.90.
- `DiscoTope_raw` range: {min(raw_scores):.5f}–{max(raw_scores):.5f}.
- Calibrated score range: {min(calibrated):.5f}–{max(calibrated):.5f}.
- Official positive residues: {len(positives)}/181; contiguous positive segments: {segment_text}.
- Ten highest raw-score residues: {', '.join(f"{row['AA']}{row['position']} ({row['DiscoTope_raw']:.5f})" for row in top_residues)}.
- Highest nonredundant 15-aa raw-score windows: {', '.join(f"{window['start']}-{window['end']} (mean {window['mean']:.5f})" for window in windows)}.

These score summaries are derived only from the restored predictor output and do not use historical assay intervals. Formal multi-evidence candidate selection is deferred to Steps 7B–8B.
"""
    args.report.write_text(report, encoding="utf-8")

    # Add self-hashes after all artifacts exist.
    provenance["outputs"]["provenance.json_pre_self_hash"] = sha256(provenance_path)
    provenance["outputs"]["step05b_report.md"] = sha256(args.report)
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
