#!/usr/bin/env python3
"""Phase 1B candidate reranking with restored conformational evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
from itertools import combinations
from pathlib import Path

import pandas as pd

try:
    from scripts.candidates import (
        _best_family_match,
        _find_reference_id,
        bootstrap_region_conservation,
        cluster_windows,
        overlap_coefficient,
        sliding_windows,
    )
    from scripts.conservation import parse_fasta
    from scripts.integration import read_fasta
except ModuleNotFoundError:
    from candidates import (
        _best_family_match,
        _find_reference_id,
        bootstrap_region_conservation,
        cluster_windows,
        overlap_coefficient,
        sliding_windows,
    )
    from conservation import parse_fasta
    from integration import read_fasta


SCALES = (15, 20, 25)
LOO_SETTINGS = ("Full", "-Epitope", "-Surface", "-Conservation")
SENSITIVITY_SETTINGS = ("BepiPred-only", "DiscoTope-only")


def score_settings_phase1b(table: pd.DataFrame) -> dict[str, list[float]]:
    b = table["BepiPred_percentile"].astype(float)
    d = table["DiscoTope_percentile"].astype(float)
    e = table["EpitopeConsensus"].astype(float)
    r = table["RSA_percentile"].astype(float)
    c = table["Conservation_percentile"].astype(float)
    settings = {
        "Full": ((e + r + c) / 3.0).tolist(),
        "-Epitope": ((r + c) / 2.0).tolist(),
        "-Surface": ((e + c) / 2.0).tolist(),
        "-Conservation": ((e + r) / 2.0).tolist(),
        "BepiPred-only": ((b + r + c) / 3.0).tolist(),
        "DiscoTope-only": ((d + r + c) / 3.0).tolist(),
        "4way-sensitivity": ((b + d + r + c) / 4.0).tolist(),
    }
    if "FinalScore_Phase1B" in table:
        error = max(
            abs(a - b)
            for a, b in zip(settings["Full"], table["FinalScore_Phase1B"])
        )
        if error >= 1e-12:
            raise ValueError(f"Full setting disagrees with Phase 1B score: {error}")
    if "FinalScore_Phase1A" in table:
        error = max(
            abs(a - b)
            for a, b in zip(
                settings["BepiPred-only"], table["FinalScore_Phase1A"]
            )
        )
        if error >= 1e-12:
            raise ValueError(f"BepiPred-only disagrees with Phase 1A score: {error}")
    return settings


def rewrite_phase1b_pdb_bfactors(
    source: Path,
    target: Path,
    scores: dict[int, float],
    *,
    label: str,
) -> None:
    seen: set[int] = set()
    output = [
        f"REMARK 900 B-FACTOR STORES {label} X 100\n",
        "REMARK 900 PHASE 1B PRE-UNBLIND; DISCOTOPE AVAILABLE TRUE\n",
    ]
    for line in source.read_text(encoding="ascii").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == "A":
            position = int(line[22:26])
            if position in scores:
                seen.add(position)
                newline = "\n" if line.endswith("\n") else ""
                body = line.rstrip("\r\n")
                if len(body) < 66:
                    body = body.ljust(66)
                line = (
                    f"{body[:60]}{float(scores[position]) * 100.0:6.2f}"
                    f"{body[66:]}{newline}"
                )
        output.append(line)
    missing = sorted(set(scores) - seen)
    if missing:
        raise ValueError(f"PDB is missing project residues: {missing}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("".join(output), encoding="ascii")


def _family_collections(settings: dict[str, list[float]]):
    windows_by_setting: dict[str, list[dict]] = {}
    families_by_setting: dict[str, list[dict]] = {}
    full_by_scale: dict[int, list[dict]] = {}
    for setting, scores in settings.items():
        combined: list[dict] = []
        for scale in SCALES:
            windows = sliding_windows(scores, scale, setting)
            ordered = sorted(
                windows, key=lambda row: (-row["mean_score"], row["start"])
            )
            for rank, window in enumerate(ordered, start=1):
                window["rank_within_scale"] = rank
            combined.extend(windows)
            if setting == "Full":
                full_by_scale[scale] = windows
        windows_by_setting[setting] = combined
        families_by_setting[setting] = cluster_windows(combined)
    return windows_by_setting, families_by_setting, full_by_scale


def _stability_matrix(
    top_families: list[dict], families_by_setting: dict[str, list[dict]]
) -> pd.DataFrame:
    rows: list[dict] = []
    for candidate_rank, family in enumerate(top_families, start=1):
        candidate_id = f"C{candidate_rank}-B"
        candidate = family["representative"]
        for setting in LOO_SETTINGS:
            match = _best_family_match(candidate, families_by_setting[setting])
            representative = match["representative"]
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "setting": setting,
                    "ranking": match["rank"],
                    "mean_score": representative["mean_score"],
                    "region_overlap": match["overlap"],
                    "matched_start": representative["start"],
                    "matched_end": representative["end"],
                    "matched_length": representative["length"],
                    "matched_scale": representative["scale"],
                    "stable_top5_overlap_gt_0_5": (
                        match["rank"] <= 5 and match["overlap"] > 0.5
                    ),
                }
            )
    return pd.DataFrame(rows)


def _scale_matches(
    top_families: list[dict], full_by_scale: dict[int, list[dict]]
) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for candidate_rank, family in enumerate(top_families, start=1):
        candidate_id = f"C{candidate_rank}-B"
        candidate = family["representative"]
        stable = []
        for scale in SCALES:
            scale_families = cluster_windows(full_by_scale[scale])[:5]
            match = _best_family_match(candidate, scale_families)
            if match["rank"] <= 5 and match["overlap"] > 0.5:
                stable.append(scale)
        result[candidate_id] = stable
    return result


def _candidate_table(
    table: pd.DataFrame,
    top_families: list[dict],
    stability: pd.DataFrame,
    scale_matches: dict[str, list[int]],
    four_way: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for rank, family in enumerate(top_families, start=1):
        representative = family["representative"]
        candidate_id = f"C{rank}-B"
        segment = table.iloc[representative["start"] - 1 : representative["end"]]
        stable_settings = int(
            stability[
                (stability["candidate_id"] == candidate_id)
                & stability["stable_top5_overlap_gt_0_5"]
            ].shape[0]
        )
        scales = scale_matches[candidate_id]
        fw = four_way[four_way["candidate_id"] == candidate_id].iloc[0]
        median_plddt = float(segment["pLDDT"].median())
        rows.append(
            {
                "candidate_id": candidate_id,
                "family_rank": rank,
                "start": representative["start"],
                "end": representative["end"],
                "length": representative["length"],
                "representative_scale": representative["scale"],
                "FinalScore_Phase1B": float(
                    segment["FinalScore_Phase1B"].mean()
                ),
                "FinalScore_4way_sensitivity_mean": float(
                    segment["FinalScore_4way_sensitivity"].mean()
                ),
                "BepiPred_mean": float(segment["BepiPred_raw"].mean()),
                "BepiPred_percentile_mean": float(
                    segment["BepiPred_percentile"].mean()
                ),
                "DiscoTope_mean": float(segment["DiscoTope_raw"].mean()),
                "DiscoTope_percentile_mean": float(
                    segment["DiscoTope_percentile"].mean()
                ),
                "EpitopeConsensus_mean": float(
                    segment["EpitopeConsensus"].mean()
                ),
                "RSA_mean": float(segment["RSA_raw"].mean()),
                "RSA_percentile_mean": float(
                    segment["RSA_percentile"].mean()
                ),
                "Conservation_mean": float(
                    segment["Conservation_raw"].mean()
                ),
                "Conservation_percentile_mean": float(
                    segment["Conservation_percentile"].mean()
                ),
                "median_pLDDT": median_plddt,
                "fraction_pLDDT_ge_70": float((segment["pLDDT"] >= 70).mean()),
                "structural_confidence": "HIGH" if median_plddt >= 70 else "LOW",
                "scale_stability": f"{len(scales)}/3",
                "scale_stability_fraction": len(scales) / 3.0,
                "stable_scales": ";".join(str(value) for value in scales),
                "LOO_stability": f"{stable_settings}/4",
                "LOO_stability_fraction": stable_settings / 4.0,
                "four_way_sensitivity_rank": int(fw["ranking"]),
                "four_way_sensitivity_overlap": float(fw["region_overlap"]),
                "family_member_windows": len(family["members"]),
            }
        )
    return pd.DataFrame(rows)


def _matching_table(
    top_families: list[dict],
    target_families: list[dict],
    *,
    setting: str,
) -> pd.DataFrame:
    rows = []
    for rank, family in enumerate(top_families, start=1):
        candidate_id = f"C{rank}-B"
        candidate = family["representative"]
        match = _best_family_match(candidate, target_families)
        representative = match["representative"]
        rows.append(
            {
                "candidate_id": candidate_id,
                "setting": setting,
                "ranking": match["rank"],
                "mean_score": representative["mean_score"],
                "region_overlap": match["overlap"],
                "matched_start": representative["start"],
                "matched_end": representative["end"],
                "matched_scale": representative["scale"],
            }
        )
    return pd.DataFrame(rows)


def _window_frame(table: pd.DataFrame, windows: list[dict]) -> pd.DataFrame:
    rows = []
    for window in windows:
        segment = table.iloc[window["start"] - 1 : window["end"]]
        rows.append(
            {
                **window,
                "BepiPred_percentile_mean": segment[
                    "BepiPred_percentile"
                ].mean(),
                "DiscoTope_percentile_mean": segment[
                    "DiscoTope_percentile"
                ].mean(),
                "EpitopeConsensus_mean": segment["EpitopeConsensus"].mean(),
                "RSA_percentile_mean": segment["RSA_percentile"].mean(),
                "Conservation_percentile_mean": segment[
                    "Conservation_percentile"
                ].mean(),
                "FinalScore_Phase1B_mean": segment[
                    "FinalScore_Phase1B"
                ].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values("rank_within_scale")


def _epitope_sensitivity(
    table: pd.DataFrame, families_by_setting: dict[str, list[dict]], top_families
) -> pd.DataFrame:
    primary = [family["representative"] for family in top_families]
    rows = []
    for setting in ("Full", *SENSITIVITY_SETTINGS):
        for rank, family in enumerate(families_by_setting[setting][:5], start=1):
            representative = family["representative"]
            matches = [
                overlap_coefficient(
                    representative["start"],
                    representative["end"],
                    candidate["start"],
                    candidate["end"],
                )
                for candidate in primary
            ]
            best = max(range(len(matches)), key=matches.__getitem__)
            matched_primary = (
                f"C{best + 1}-B" if matches[best] > 0 else "NO_TOP5_MATCH"
            )
            segment = table.iloc[
                representative["start"] - 1 : representative["end"]
            ]
            rows.append(
                {
                    "model_setting": setting,
                    "rank": rank,
                    "start": representative["start"],
                    "end": representative["end"],
                    "length": representative["length"],
                    "scale": representative["scale"],
                    "mean_score": representative["mean_score"],
                    "matched_primary_candidate": matched_primary,
                    "overlap_with_primary": matches[best],
                    "BepiPred_percentile_mean": segment[
                        "BepiPred_percentile"
                    ].mean(),
                    "DiscoTope_percentile_mean": segment[
                        "DiscoTope_percentile"
                    ].mean(),
                }
            )
    return pd.DataFrame(rows)


def _phase_comparison(
    phase1a: pd.DataFrame,
    phase1b: pd.DataFrame,
    feature_table: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for a in phase1a.itertuples(index=False):
        overlaps = [
            overlap_coefficient(a.start, a.end, int(b.start), int(b.end))
            for b in phase1b.itertuples(index=False)
        ]
        best_index = max(range(len(overlaps)), key=overlaps.__getitem__)
        b = list(phase1b.itertuples(index=False))[best_index]
        segment_a = feature_table.iloc[int(a.start) - 1 : int(a.end)]
        segment_b = feature_table.iloc[int(b.start) - 1 : int(b.end)]
        if overlaps[best_index] <= 0:
            rows.append(
                {
                    "phase1a_candidate": f"{a.candidate_id}-A",
                    "phase1a_start": int(a.start),
                    "phase1a_end": int(a.end),
                    "rank_A": int(a.family_rank),
                    "phase1a_primary_score": float(a.FinalScore),
                    "phase1b_candidate": "NO_TOP5_MATCH",
                    "phase1b_start": pd.NA,
                    "phase1b_end": pd.NA,
                    "rank_B": pd.NA,
                    "phase1b_primary_score": math.nan,
                    "region_overlap": 0.0,
                    "rank_change_A_minus_B": pd.NA,
                    "mean_3d_contribution_delta_A_region": float(
                        (
                            segment_a["FinalScore_Phase1B"]
                            - segment_a["FinalScore_Phase1A"]
                        ).mean()
                    ),
                    "mean_3d_contribution_delta_B_region": math.nan,
                    "DiscoTope_percentile_mean_A_region": float(
                        segment_a["DiscoTope_percentile"].mean()
                    ),
                    "DiscoTope_percentile_mean_B_region": math.nan,
                }
            )
            continue
        rows.append(
            {
                "phase1a_candidate": f"{a.candidate_id}-A",
                "phase1a_start": int(a.start),
                "phase1a_end": int(a.end),
                "rank_A": int(a.family_rank),
                "phase1a_primary_score": float(a.FinalScore),
                "phase1b_candidate": b.candidate_id,
                "phase1b_start": int(b.start),
                "phase1b_end": int(b.end),
                "rank_B": int(b.family_rank),
                "phase1b_primary_score": float(b.FinalScore_Phase1B),
                "region_overlap": overlaps[best_index],
                "rank_change_A_minus_B": int(a.family_rank) - int(b.family_rank),
                "mean_3d_contribution_delta_A_region": float(
                    (
                        segment_a["FinalScore_Phase1B"]
                        - segment_a["FinalScore_Phase1A"]
                    ).mean()
                ),
                "mean_3d_contribution_delta_B_region": float(
                    (
                        segment_b["FinalScore_Phase1B"]
                        - segment_b["FinalScore_Phase1A"]
                    ).mean()
                ),
                "DiscoTope_percentile_mean_A_region": float(
                    segment_a["DiscoTope_percentile"].mean()
                ),
                "DiscoTope_percentile_mean_B_region": float(
                    segment_b["DiscoTope_percentile"].mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def _plot_candidate_profile(path: Path, table: pd.DataFrame, candidates: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tracks = [
        ("BepiPred_percentile", "#1565C0", "BepiPred"),
        ("DiscoTope_percentile", "#6A1B9A", "DiscoTope"),
        ("RSA_percentile", "#00897B", "RSA"),
        ("Conservation_percentile", "#EF6C00", "Conservation"),
        ("FinalScore_Phase1B", "#C62828", "FinalScore Phase 1B"),
    ]
    colors = ["#D32F2F", "#1976D2", "#388E3C", "#7B1FA2", "#F57C00"]
    fig, axes = plt.subplots(5, 1, figsize=(13, 10.5), sharex=True)
    for axis, (column, color, label) in zip(axes, tracks):
        axis.plot(table["position"], table[column], color=color, linewidth=1.3)
        axis.set_ylabel(label)
        axis.set_ylim(-0.03, 1.03)
        axis.grid(axis="y", alpha=0.18)
        for span_color, row in zip(colors, candidates.itertuples(index=False)):
            axis.axvspan(row.start, row.end, color=span_color, alpha=0.08)
        axis.spines[["top", "right"]].set_visible(False)
    for span_color, row in zip(colors, candidates.itertuples(index=False)):
        axes[-1].axvspan(
            row.start,
            row.end,
            color=span_color,
            alpha=0.18,
            label=f"{row.candidate_id}: {row.start}-{row.end}",
        )
    axes[-1].legend(frameon=False, ncol=3, fontsize=8, loc="lower left")
    axes[-1].set_xlabel("Project position (1-181)")
    axes[-1].set_xlim(1, 181)
    fig.suptitle("IL-24 Phase 1B evidence tracks and primary candidate families")
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_robustness(path: Path, stability: pd.DataFrame, candidates: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matrix = []
    overlaps = []
    for candidate_id in candidates["candidate_id"]:
        subset = stability[stability["candidate_id"] == candidate_id].set_index(
            "setting"
        )
        matrix.append([subset.loc[value, "ranking"] for value in LOO_SETTINGS])
        overlaps.append(
            [subset.loc[value, "region_overlap"] for value in LOO_SETTINGS]
        )
    values = np.asarray(matrix, dtype=float)
    fig, axis = plt.subplots(figsize=(8.8, 5.5))
    image = axis.imshow(
        values,
        cmap="viridis_r",
        aspect="auto",
        vmin=1,
        vmax=max(5, values.max()),
    )
    axis.set_xticks(range(4), labels=LOO_SETTINGS)
    axis.set_yticks(range(len(candidates)), labels=candidates["candidate_id"])
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            axis.text(
                column,
                row,
                f"rank {int(values[row, column])}\nov {overlaps[row][column]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )
    axis.set_title("Phase 1B grouped leave-one-evidence-out robustness")
    fig.colorbar(image, ax=axis, label="Matched family rank")
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_rank_change(path: Path, comparison: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8.5, 6.2))
    for index, row in comparison.iterrows():
        unmatched = pd.isna(row["rank_B"])
        rank_b = 6 if unmatched else row["rank_B"]
        axis.plot([0, 1], [row["rank_A"], rank_b], marker="o", linewidth=1.6)
        axis.text(-0.03, row["rank_A"], row["phase1a_candidate"], ha="right", va="center")
        axis.text(1.03, rank_b, row["phase1b_candidate"], ha="left", va="center")
    axis.set_xticks([0, 1], labels=["Phase 1A", "Phase 1B matched family"])
    axis.set_yticks(range(1, 7), labels=["1", "2", "3", "4", "5", "not Top 5"])
    axis.set_ylim(6.5, 0.5)
    axis.set_xlim(-0.3, 1.3)
    axis.set_ylabel("Candidate rank")
    axis.set_title("Phase 1A to Phase 1B candidate-family rank mapping")
    axis.grid(axis="y", alpha=0.2)
    axis.spines[["top", "right", "bottom"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_phase1a_manifest(root: Path, manifest: Path) -> dict:
    rows = list(csv.DictReader(manifest.open(encoding="utf-8-sig", newline="")))
    mismatches = []
    for row in rows:
        relative = row["path"].replace("\\", "/")
        path = root / relative
        observed = _sha256(path) if path.is_file() else None
        if observed != row["sha256"]:
            mismatches.append(
                {"path": relative, "expected": row["sha256"], "observed": observed}
            )
    return {
        "manifest": str(manifest),
        "manifest_sha256": _sha256(manifest),
        "entries": len(rows),
        "mismatches": mismatches,
        "all_unchanged": not mismatches,
    }


def _super_family_notes(candidates: pd.DataFrame) -> list[str]:
    notes = []
    for left, right in combinations(candidates.itertuples(index=False), 2):
        intersection = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
        gap = max(left.start, right.start) - min(left.end, right.end) - 1
        if intersection > 0 or 0 <= gap <= 5:
            notes.append(
                f"{left.candidate_id} ({left.start}-{left.end}) and "
                f"{right.candidate_id} ({right.start}-{right.end}) may be interpreted "
                "as a broader biological hotspot super-family; formal families remain separate."
            )
    return notes


def _write_report(
    path: Path,
    candidates: pd.DataFrame,
    stability: pd.DataFrame,
    sensitivity: pd.DataFrame,
    comparison: pd.DataFrame,
    correlations: dict,
    audit: dict,
    super_notes: list[str],
    discotope_provenance: dict,
):
    candidate_lines = [
        "| Candidate | Region | FinalB | BepiPred | DiscoTope | Epitope consensus | RSA | Conservation | median pLDDT | pLDDT>=70 | Scale | LOO | Bootstrap median [95% CI] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in candidates.itertuples(index=False):
        candidate_lines.append(
            f"| {row.candidate_id} | {row.start}-{row.end} | {row.FinalScore_Phase1B:.4f} | "
            f"{row.BepiPred_mean:.4f} | {row.DiscoTope_mean:.4f} | {row.EpitopeConsensus_mean:.4f} | "
            f"{row.RSA_mean:.4f} | {row.Conservation_mean:.4f} | {row.median_pLDDT:.2f} | "
            f"{row.fraction_pLDDT_ge_70:.2f} | {row.scale_stability} ({row.stable_scales}) | "
            f"{row.LOO_stability} | {row.conservation_bootstrap_median:.4f} "
            f"[{row.conservation_bootstrap_ci_lower:.4f}, {row.conservation_bootstrap_ci_upper:.4f}] |"
        )
    stability_lines = [
        "| Candidate | Full | -Epitope | -Surface | -Conservation |",
        "|---|---:|---:|---:|---:|",
    ]
    for candidate_id in candidates["candidate_id"]:
        subset = stability[stability["candidate_id"] == candidate_id].set_index("setting")
        cells = []
        for setting in LOO_SETTINGS:
            row = subset.loc[setting]
            cells.append(
                f"rank {int(row['ranking'])}; mean {row['mean_score']:.4f}; ov {row['region_overlap']:.2f}"
            )
        stability_lines.append(f"| {candidate_id} | " + " | ".join(cells) + " |")
    compare_lines = [
        "| Phase 1A | Phase 1B match | Overlap | Rank A -> B | 3D delta (A region) |",
        "|---|---|---:|---:|---:|",
    ]
    for row in comparison.itertuples(index=False):
        if pd.isna(row.rank_B):
            compare_lines.append(
                f"| {row.phase1a_candidate} {row.phase1a_start}-{row.phase1a_end} | "
                f"NO_TOP5_MATCH | 0.00 | {row.rank_A} -> not Top 5 | "
                f"{row.mean_3d_contribution_delta_A_region:+.4f} |"
            )
        else:
            compare_lines.append(
                f"| {row.phase1a_candidate} {row.phase1a_start}-{row.phase1a_end} | "
                f"{row.phase1b_candidate} {int(row.phase1b_start)}-{int(row.phase1b_end)} | "
                f"{row.region_overlap:.2f} | {row.rank_A} -> {int(row.rank_B)} | "
                f"{row.mean_3d_contribution_delta_A_region:+.4f} |"
            )
    sensitivity_summary = []
    for setting in SENSITIVITY_SETTINGS:
        subset = sensitivity[sensitivity["model_setting"] == setting]
        sensitivity_summary.append(
            f"- {setting} Top 5: "
            + ", ".join(
                f"{int(row.start)}-{int(row.end)} (rank {int(row['rank'])}, "
                f"primary match {row.matched_primary_candidate}, ov {row.overlap_with_primary:.2f})"
                for _, row in subset.iterrows()
            )
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""# IL-24 Phase 1B Pre-unblinding Report

## 1. Why Phase 1B was performed

The preregistered DiscoTope-3.0 structural-epitope branch was missing in Phase 1A because the allowed installation attempts ended in external dependency failure. Before any historical experimental unblinding, the user explicitly authorized a post-freeze restoration of that missing branch. This is a protocol amendment, not a retroactive claim that the restoration was part of the original successful run.

## 2. Frozen Phase 1A

Phase 1A outputs remained read-only. The original freeze manifest was replayed after Phase 1B: {audit['entries']}/{audit['entries']} entries matched, zero mismatches, manifest SHA-256 `{audit['manifest_sha256']}`. Historical experimental regions were not accessed.

## 3. DiscoTope-3.0

Official DiscoTope-3.0 commit `35d9f2e55f97eaba2a7acefbc394db58fb9670bc` ran in `il24-dt3-phase1b` (Python 3.14.7, Torch 2.10.0/CUDA 12.8). B2 added only GCC/G++ to compile the official `biotraj` dependency; predictor code, weights, and requirement versions were unchanged. Official ESM-IF1, 100 XGBoost models, and two GAM calibrators loaded successfully. AlphaFold mode used the RTX 4050 GPU without OOM or CPU fallback; wall time was {discotope_provenance['inference']['wall_seconds']:.2f} s.

DiscoTope is a structure-aware B-cell conformational epitope predictor. Unlike a sequence-only model, it encodes the folded backbone so residues far apart in the linear chain can jointly contribute to an antibody-accessible 3D surface patch. RSA asks whether a residue is exposed; DiscoTope asks whether that residue and its 3D neighborhood resemble an epitope.

## 4. 3D epitope results

The official output mapped exactly to positions 1-181 with no non-finite values. `DiscoTope_raw` is the official uncalibrated continuous score used for percentile ranking; official calibrated scores and the default 0.90 classification were separately preserved. Raw range: 0.00235-0.57638; official positive calls: 38/181.

## 5. BepiPred vs DiscoTope

Raw Pearson r = {correlations['pearson_raw']:.4f}; raw Spearman rho = {correlations['spearman_raw']:.4f}. The moderate agreement indicates both convergence and complementary information rather than redundancy.

## 6. Full evidence integration

All features use the Phase 1A average-tie percentile transform on [0,1]. `E=(B+D)/2`, then the primary score is `FinalScore_Phase1B=(E+R+C)/3`, giving one-third weight each to epitope consensus, surface accessibility, and evolution. `(B+D+R+C)/4` is sensitivity-only. pLDDT is a confidence flag and DSSP is annotation; neither changes the score.

## 7. Phase 1B candidate regions

All 15-, 20-, and 25-aa windows were freshly recomputed from FinalScore Phase 1B. The unchanged family rule is overlap coefficient `intersection / shorter length > 0.5`; no single scale or unique Top1 was selected.

{chr(10).join(candidate_lines)}

## 8. Robustness

Grouped leave-one-evidence-out results are reported as matched family rank, matched mean, and overlap. Scale stability uses region-level matching across all three fixed scales. Conservation uncertainty reuses the frozen Step 6 alignment with exactly 100 bootstrap repetitions.

{chr(10).join(stability_lines)}

Epitope-model sensitivity:

{chr(10).join(sensitivity_summary)}

The BepiPred-only residue score reproduces Phase 1A to numerical precision. The DiscoTope-only and combined rankings were computed independently; no weights were tuned to preserve a Phase 1A hotspot.

## 9. Phase 1A vs Phase 1B

{chr(10).join(compare_lines)}

## 10. Interpretation

Regions with high mean BepiPred and DiscoTope percentile have convergent sequence-based and structure-based epitope support. Regions reduced by Phase 1B had strong sequence evidence but weaker conformational support; newly promoted regions gained structure-supported priority. These are all acceptable pre-unblinding outcomes.

{chr(10).join('- ' + note for note in super_notes) if super_notes else '- No adjacent formally split Top-5 families required super-family interpretation.'}

Super-family language is interpretive only: it does not merge candidates or alter the fixed clustering threshold.

## 11. Limitations

- These are not experimentally validated epitopes.
- These are not receptor interfaces.
- These are not confirmed neutralizing epitopes.
- These are not periodontitis-specific functional interfaces.
- AlphaFold confidence and computational predictor scores do not establish biological function.

## 12. Ready for unblinding

Phase 1B pre-unblinding analysis frozen. Historical experimental intervals remain blinded; Step 9 was not performed.
""",
        encoding="utf-8",
    )


def run(args: argparse.Namespace) -> None:
    table = pd.read_csv(args.feature_table)
    reference = read_fasta(args.reference)
    if len(table) != 181 or table["position"].tolist() != list(range(1, 182)):
        raise ValueError("Phase 1B feature table must contain positions 1-181")
    if "".join(table["AA"].astype(str)) != reference:
        raise ValueError("Phase 1B feature table does not match frozen reference")
    required = [
        "pLDDT",
        "BepiPred_raw",
        "BepiPred_percentile",
        "DiscoTope_raw",
        "DiscoTope_percentile",
        "EpitopeConsensus",
        "RSA_raw",
        "RSA_percentile",
        "Conservation_raw",
        "Conservation_percentile",
        "FinalScore_Phase1A",
        "FinalScore_Phase1B",
        "FinalScore_4way_sensitivity",
    ]
    if not table[required].map(lambda value: math.isfinite(float(value))).all().all():
        raise ValueError("Phase 1B feature table contains NaN or infinity")

    settings = score_settings_phase1b(table)
    windows_by_setting, families_by_setting, full_by_scale = _family_collections(settings)
    top_families = families_by_setting["Full"][:5]
    stability = _stability_matrix(top_families, families_by_setting)
    scales = _scale_matches(top_families, full_by_scale)
    four_way = _matching_table(
        top_families,
        families_by_setting["4way-sensitivity"],
        setting="4way-sensitivity",
    )
    candidates = _candidate_table(table, top_families, stability, scales, four_way)

    alignment = parse_fasta(args.alignment)
    reference_id = _find_reference_id(alignment, reference)
    bootstrap = pd.DataFrame(
        bootstrap_region_conservation(
            alignment,
            reference_id,
            candidates[["candidate_id", "start", "end"]].to_dict("records"),
            replicates=100,
            seed=2406,
        )
    )
    candidates = candidates.merge(bootstrap, on="candidate_id", validate="one_to_one")
    sensitivity = _epitope_sensitivity(table, families_by_setting, top_families)
    phase1a = pd.read_csv(args.phase1a_candidates)
    comparison = _phase_comparison(phase1a, candidates, table)
    audit = audit_phase1a_manifest(args.root, args.phase1a_manifest)
    if not audit["all_unchanged"]:
        raise RuntimeError(f"Frozen Phase 1A hash mismatch: {audit['mismatches']}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scale in SCALES:
        _window_frame(table, full_by_scale[scale]).to_csv(
            args.output_dir / f"window_{scale}_scores_3d.csv", index=False
        )
    candidates.to_csv(
        args.output_dir / "candidate_regions_phase1b_preunblind.csv", index=False
    )
    stability.to_csv(
        args.output_dir / "candidate_stability_matrix_phase1b.csv", index=False
    )
    sensitivity.to_csv(args.output_dir / "epitope_model_sensitivity.csv", index=False)
    four_way.to_csv(
        args.output_dir / "four_way_weighting_sensitivity.csv", index=False
    )
    comparison.to_csv(
        args.output_dir / "phase1a_vs_phase1b_candidates.csv", index=False
    )
    (args.output_dir / "phase1a_hash_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    _plot_candidate_profile(
        args.output_dir / "02_candidate_profile_phase1b.png", table, candidates
    )
    _plot_robustness(
        args.output_dir / "03_robustness_heatmap_phase1b.png",
        stability,
        candidates,
    )
    _plot_rank_change(
        args.output_dir / "06_phase1a_vs_phase1b_candidate_rank.png", comparison
    )
    shutil.copyfile(
        args.integration_dir / "01_all_evidence_tracks_phase1b.png",
        args.output_dir / "01_all_evidence_tracks_phase1b.png",
    )
    shutil.copyfile(
        args.integration_dir / "04_bepipred_vs_discotope.png",
        args.output_dir / "04_bepipred_vs_discotope.png",
    )
    shutil.copyfile(
        args.integration_dir / "05_phase1a_vs_phase1b_profile.png",
        args.output_dir / "05_phase1a_vs_phase1b_profile.png",
    )
    rewrite_phase1b_pdb_bfactors(
        args.pdb,
        args.output_dir / "IL24_181_phase1b_candidate_score.pdb",
        dict(zip(table["position"].astype(int), table["FinalScore_Phase1B"])),
        label="PHASE1B FINALSCORE",
    )
    rewrite_phase1b_pdb_bfactors(
        args.pdb,
        args.output_dir / "IL24_181_discotope_score.pdb",
        dict(zip(table["position"].astype(int), table["DiscoTope_percentile"])),
        label="DISCOTOPE PERCENTILE",
    )

    correlations = json.loads(args.correlations.read_text(encoding="utf-8"))
    discotope_provenance = json.loads(
        args.discotope_provenance.read_text(encoding="utf-8")
    )
    super_notes = _super_family_notes(candidates)
    _write_report(
        args.report,
        candidates,
        stability,
        sensitivity,
        comparison,
        correlations,
        audit,
        super_notes,
        discotope_provenance,
    )
    (args.output_dir / "provenance.json").write_text(
        json.dumps(
            {
                "phase": "Phase 1B pre-unblind",
                "historical_region_still_blinded": True,
                "primary_score": "(((B+D)/2)+R+C)/3",
                "window_scales": list(SCALES),
                "clustering": "greedy score order; overlap coefficient > 0.5",
                "top_family_count": 5,
                "LOO_settings": list(LOO_SETTINGS),
                "sensitivity_settings": [*SENSITIVITY_SETTINGS, "4way-sensitivity"],
                "bootstrap_repetitions": 100,
                "bootstrap_seed": 2406,
                "phase1a_hashes_unchanged": True,
                "step9_performed": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--feature-table", type=Path, default=Path("results/07b_integration_with_3d/residue_feature_table_3d.csv"))
    parser.add_argument("--integration-dir", type=Path, default=Path("results/07b_integration_with_3d"))
    parser.add_argument("--reference", type=Path, default=Path("inputs/reference_181.fasta"))
    parser.add_argument("--alignment", type=Path, default=Path("results/06_conservation/mafft_alignment.fasta"))
    parser.add_argument("--pdb", type=Path, default=Path("inputs/AF-Q925S4-F1-model_v6.pdb"))
    parser.add_argument("--phase1a-candidates", type=Path, default=Path("results/phase1a_baseline/candidate_regions_preunblind.csv"))
    parser.add_argument("--phase1a-manifest", type=Path, default=Path("logs/phase1_preunblind_freeze_sha256.csv"))
    parser.add_argument("--correlations", type=Path, default=Path("results/07b_integration_with_3d/bepipred_discotope_correlations.json"))
    parser.add_argument("--discotope-provenance", type=Path, default=Path("results/05b_discotope_restored/provenance.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/08b_candidates_with_3d"))
    parser.add_argument("--report", type=Path, default=Path("reports/phase1b_preunblind_report.md"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
