#!/usr/bin/env python
"""Extract robust IL-24 candidate region families from the frozen Step 7 table."""

from __future__ import annotations

import argparse
import math
import random
import statistics
from pathlib import Path

import pandas as pd

if __package__:
    from scripts.conservation import (
        _column_conservation,
        _percentile,
        _reference_columns,
        parse_fasta,
    )
    from scripts.integration import read_fasta
else:
    from conservation import _column_conservation, _percentile, _reference_columns, parse_fasta
    from integration import read_fasta


SCALES = (15, 20, 25)
SETTINGS = ("Full", "-Epitope", "-Surface", "-Conservation")


def overlap_coefficient(start_a: int, end_a: int, start_b: int, end_b: int) -> float:
    intersection = max(0, min(end_a, end_b) - max(start_a, start_b) + 1)
    shorter = min(end_a - start_a + 1, end_b - start_b + 1)
    return intersection / shorter


def sliding_windows(scores: list[float], window_size: int, setting: str) -> list[dict]:
    if window_size < 1 or window_size > len(scores):
        raise ValueError(f"Invalid window size {window_size} for {len(scores)} scores")
    windows = []
    for offset in range(len(scores) - window_size + 1):
        segment = scores[offset : offset + window_size]
        windows.append(
            {
                "start": offset + 1,
                "end": offset + window_size,
                "length": window_size,
                "scale": window_size,
                "setting": setting,
                "mean_score": statistics.fmean(segment),
            }
        )
    return windows


def cluster_windows(windows: list[dict]) -> list[dict]:
    """Greedy score-ordered family assignment against fixed representatives."""
    ordered = sorted(
        windows,
        key=lambda row: (-row["mean_score"], row["start"], row["end"], row["scale"]),
    )
    families: list[dict] = []
    for window in ordered:
        for family in families:
            representative = family["representative"]
            if (
                overlap_coefficient(
                    window["start"],
                    window["end"],
                    representative["start"],
                    representative["end"],
                )
                > 0.5
            ):
                family["members"].append(window)
                break
        else:
            families.append({"representative": window, "members": [window]})
    return families


def score_settings(table: pd.DataFrame) -> dict[str, list[float]]:
    b = table["BepiPred_percentile"].astype(float)
    r = table["RSA_percentile"].astype(float)
    c = table["Conservation_percentile"].astype(float)
    return {
        "Full": ((b + r + c) / 3.0).tolist(),
        "-Epitope": ((r + c) / 2.0).tolist(),
        "-Surface": ((b + c) / 2.0).tolist(),
        "-Conservation": ((b + r) / 2.0).tolist(),
    }


def _best_family_match(candidate: dict, families: list[dict]) -> dict:
    matches = []
    for rank, family in enumerate(families, start=1):
        representative = family["representative"]
        overlap = overlap_coefficient(
            candidate["start"], candidate["end"], representative["start"], representative["end"]
        )
        matches.append((overlap, -rank, rank, representative))
    overlap, _, rank, representative = max(matches)
    return {"overlap": overlap, "rank": rank, "representative": representative}


def bootstrap_region_conservation(
    alignment: dict[str, str],
    reference_id: str,
    regions: list[dict],
    replicates: int = 100,
    seed: int = 2406,
) -> list[dict]:
    if replicates != 100:
        raise ValueError("The frozen protocol requires exactly 100 bootstrap repetitions")
    orthologs = [sequence for key, sequence in alignment.items() if key != reference_id]
    if not orthologs:
        raise ValueError("Alignment contains no orthologs")
    columns = _reference_columns(alignment, reference_id)
    rng = random.Random(seed)
    region_values = {region["candidate_id"]: [] for region in regions}
    for _ in range(replicates):
        sampled = [orthologs[rng.randrange(len(orthologs))] for _ in orthologs]
        residue_scores = [
            _column_conservation([sequence[column] for sequence in sampled])[0]
            for column in columns
        ]
        for region in regions:
            values = residue_scores[region["start"] - 1 : region["end"]]
            region_values[region["candidate_id"]].append(statistics.median(values))
    results = []
    for region in regions:
        values = region_values[region["candidate_id"]]
        results.append(
            {
                "candidate_id": region["candidate_id"],
                "conservation_bootstrap_median": statistics.median(values),
                "conservation_bootstrap_ci_lower": _percentile(values, 0.025),
                "conservation_bootstrap_ci_upper": _percentile(values, 0.975),
                "bootstrap_repetitions": replicates,
                "bootstrap_seed": seed,
            }
        )
    return results


def rewrite_pdb_bfactors(source: Path, target: Path, scores: dict[int, float]) -> None:
    seen: set[int] = set()
    output = [
        "REMARK 900 B-FACTOR STORES AVAILABLE-EVIDENCE FINALSCORE X 100\n",
        "REMARK 900 DISCOTOPE_AVAILABLE FALSE; NO SUBSTITUTE MODEL USED\n",
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
                line = f"{body[:60]}{scores[position] * 100.0:6.2f}{body[66:]}{newline}"
        output.append(line)
    missing = sorted(set(scores) - seen)
    if missing:
        raise ValueError(f"PDB is missing project residues: {missing}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("".join(output), encoding="ascii")


def _find_reference_id(alignment: dict[str, str], reference: str) -> str:
    matches = [key for key, sequence in alignment.items() if sequence.replace("-", "") == reference]
    if len(matches) != 1:
        raise ValueError(f"Expected one exact project reference in alignment; found {matches}")
    return matches[0]


def _candidate_table(
    feature_table: pd.DataFrame,
    top_families: list[dict],
    stability: pd.DataFrame,
    scale_matches: dict[str, list[int]],
) -> pd.DataFrame:
    rows = []
    for rank, family in enumerate(top_families, start=1):
        representative = family["representative"]
        candidate_id = f"C{rank}"
        segment = feature_table.iloc[representative["start"] - 1 : representative["end"]]
        setting_rows = stability[stability["candidate_id"] == candidate_id]
        stable_settings = int(
            ((setting_rows["region_overlap"] > 0.5) & (setting_rows["ranking"] <= 5)).sum()
        )
        scales = scale_matches[candidate_id]
        median_plddt = float(segment["pLDDT"].median())
        rows.append(
            {
                "candidate_id": candidate_id,
                "family_rank": rank,
                "start": representative["start"],
                "end": representative["end"],
                "length": representative["length"],
                "representative_scale": representative["scale"],
                "FinalScore": float(segment["FinalScore"].mean()),
                "BepiPred": float(segment["BepiPred_raw"].mean()),
                "BepiPred_percentile": float(segment["BepiPred_percentile"].mean()),
                "RSA": float(segment["RSA_raw"].mean()),
                "RSA_percentile": float(segment["RSA_percentile"].mean()),
                "Conservation": float(segment["Conservation_raw"].mean()),
                "Conservation_percentile": float(segment["Conservation_percentile"].mean()),
                "median_pLDDT": median_plddt,
                "fraction_pLDDT_ge_70": float((segment["pLDDT"] >= 70.0).mean()),
                "structural_confidence": "HIGH" if median_plddt >= 70.0 else "LOW",
                "scale_stability": f"{len(scales)}/3",
                "scale_stability_fraction": len(scales) / 3.0,
                "stable_scales": ";".join(str(scale) for scale in scales),
                "leave_one_evidence_out_stability": f"{stable_settings}/4",
                "leave_one_evidence_out_stability_fraction": stable_settings / 4.0,
                "family_member_windows": len(family["members"]),
            }
        )
    return pd.DataFrame(rows)


def _plot_candidate_profile(path: Path, table: pd.DataFrame, candidates: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(12, 5.2))
    axis.plot(table["position"], table["FinalScore"], color="#3f007d", linewidth=1.4)
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    for color, row in zip(colors, candidates.itertuples(index=False)):
        axis.axvspan(row.start, row.end, color=color, alpha=0.18, label=f"{row.candidate_id}: {row.start}–{row.end}")
    axis.set(
        xlabel="Project position (1–181)",
        ylabel="FinalScore",
        xlim=(1, len(table)),
        ylim=(0, 1.02),
        title="IL-24 pre-unblinding robust candidate families",
    )
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, ncol=3, fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=220)
    plt.close(figure)


def _plot_robustness(path: Path, stability: pd.DataFrame, candidates: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    matrix = []
    overlaps = []
    for candidate_id in candidates["candidate_id"]:
        subset = stability[stability["candidate_id"] == candidate_id].set_index("setting")
        matrix.append([subset.loc[setting, "ranking"] for setting in SETTINGS])
        overlaps.append([subset.loc[setting, "region_overlap"] for setting in SETTINGS])
    values = np.asarray(matrix, dtype=float)
    figure, axis = plt.subplots(figsize=(8.2, 5.2))
    image = axis.imshow(values, cmap="viridis_r", aspect="auto", vmin=1, vmax=max(5, values.max()))
    axis.set_xticks(range(len(SETTINGS)), labels=SETTINGS)
    axis.set_yticks(range(len(candidates)), labels=candidates["candidate_id"])
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            axis.text(column, row, f"rank {int(values[row, column])}\nov {overlaps[row][column]:.2f}", ha="center", va="center", fontsize=8)
    axis.set_title("Leave-one-evidence-out family robustness")
    figure.colorbar(image, ax=axis, label="Matched family rank (lower is better)")
    figure.tight_layout()
    figure.savefig(path, dpi=220)
    plt.close(figure)


def _write_final_report(
    path: Path,
    feature_table: pd.DataFrame,
    bepipred: pd.DataFrame,
    candidates: pd.DataFrame,
    stability: pd.DataFrame,
) -> None:
    header = "| Candidate | Region | Length | Final | BepiPred | RSA | Conservation | median pLDDT | pLDDT >=70 | Structural confidence | Scale | LOO | Bootstrap median [95% CI] |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|"
    candidate_lines = [header]
    for row in candidates.itertuples(index=False):
        candidate_lines.append(
            f"| {row.candidate_id} | {row.start}–{row.end} | {row.length} | {row.FinalScore:.4f} | {row.BepiPred:.4f} | {row.RSA:.4f} | {row.Conservation:.4f} | {row.median_pLDDT:.2f} | {row.fraction_pLDDT_ge_70:.2f} | {row.structural_confidence} | {row.scale_stability} ({row.stable_scales}) | {row.leave_one_evidence_out_stability} | {row.conservation_bootstrap_median:.4f} [{row.conservation_bootstrap_ci_lower:.4f}, {row.conservation_bootstrap_ci_upper:.4f}] |"
        )
    ablation_lines = [
        "Cells are `rank / mean / overlap` for the matched region family.",
        "",
        "| Candidate | Full | -Epitope | -Surface | -Conservation |",
        "|---|---:|---:|---:|---:|",
    ]
    for candidate_id in candidates["candidate_id"]:
        subset = stability[stability["candidate_id"] == candidate_id].set_index("setting")
        cells = []
        for setting in SETTINGS:
            row = subset.loc[setting]
            cells.append(
                f"{int(row['ranking'])} / {row['mean_score']:.4f} / {row['region_overlap']:.2f}"
            )
        ablation_lines.append(f"| {candidate_id} | " + " | ".join(cells) + " |")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""# IL-24 Phase 1 Pre-unblinding Report

## 1. Research question

Prioritize continuous regions on the frozen 181-aa IL-24 project sequence using only preregistered evidence available before historical experimental unblinding.

## 2. Frozen input

The coordinate system is project positions 1–181 from `inputs/reference_181.fasta`. The AlphaFold structure and UniProt mapping input retain their Step 0 manifest hashes. No historical experimental interval or excluded image was opened or used.

## 3. Sequence QC

The frozen reference, PDB Chain A, and all downstream structural residue tables map continuously across 181 project residues. The formal Step 1 output was integrity-audited but not rerun.

## 4. Structure QC

The AlphaFold Chain A structure contains 181 mapped residues. pLDDT is used only as confidence metadata and never contributes to FinalScore.

## 5. BepiPred-3.0

Official standalone BepiPred-3.0 `vt_pred` uses ESM-2 protein-language-model representations learned from amino-acid sequence context, analogous to a language model learning contextual word patterns. The formal run used RTX 4050 CUDA without CPU fallback and completed in 84.504 s. The output contains {len(bepipred)} residue probabilities; range {bepipred['bepipred_positive_probability'].min():.6f}–{bepipred['bepipred_positive_probability'].max():.6f}. Continuous probabilities, not only threshold calls, enter the consensus.

## 6. Surface accessibility

FreeSASA RSA provides independent physical surface-accessibility evidence. DSSP secondary structure is annotation only and does not add or subtract score.

## 7. DiscoTope-3.0

`discotope_available = false` and status is `MISSING_EXTERNAL_DEPENDENCY_FAILURE`. Both preregistered official installation attempts were preserved; Attempt 2 failed during WSL DNS/PyPI dependency retrieval. No third attempt, substitute model, or artificial 3D epitope score was used. This preregistered structural-epitope evidence branch is therefore missing.

## 8. Evolutionary conservation

MAFFT-aligned mammalian orthologs were scored with normalized Shannon conservation. Candidate-level uncertainty uses exactly 100 ortholog bootstrap resamples.

## 9. Consensus methodology

The primary analysis is a transparent available-evidence consensus. BepiPred probability, RSA, and conservation were independently average-tie percentile ranked to [0,1], then combined as `FinalScore = (B + R + C) / 3`. pLDDT is a structure-confidence flag; DSSP is annotation.

## 10. Candidate regions

All 15-, 20-, and 25-aa windows were retained. Score-ordered interval families use overlap coefficient `intersection / shorter length > 0.5`; no single window size or unique Top1 was selected.

{chr(10).join(candidate_lines)}

## 11. Robustness analysis

Each Full candidate was region-matched against independently recomputed Full, -Epitope, -Surface, and -Conservation families. The stability matrix contains {len(stability)} candidate-setting records with rank, mean score, and region overlap. Scale stability is defined by overlap with top-five scale-specific families rather than exact endpoints.

{chr(10).join(ablation_lines)}

## 12. Limitations

The DiscoTope structural-epitope branch is missing, BepiPred thresholds are not protein-specific calibration, AlphaFold confidence is not functional evidence, and conservation depends on the curated ortholog set. These results are multi-evidence computational candidate regions—not experimentally validated epitopes, receptor interfaces, or confirmed neutralizing epitopes.

Phase 1 remains pre-unblinded. Step 9 was not performed.
""",
        encoding="utf-8",
    )


def run(args: argparse.Namespace) -> None:
    table = pd.read_csv(args.feature_table)
    reference = read_fasta(args.reference)
    if len(table) != 181 or table["position"].tolist() != list(range(1, 182)):
        raise ValueError("Step 7 feature table must contain project positions 1–181")
    if "".join(table["AA"]) != reference:
        raise ValueError("Step 7 AA sequence does not match reference_181.fasta")
    required_numeric = [
        "pLDDT",
        "BepiPred_raw",
        "BepiPred_percentile",
        "SASA",
        "RSA_raw",
        "RSA_percentile",
        "Conservation_raw",
        "Conservation_percentile",
        "FinalScore",
    ]
    if not table[required_numeric].map(lambda value: math.isfinite(float(value))).all().all():
        raise ValueError("Step 7 feature table contains NaN or infinity")

    settings = score_settings(table)
    windows_by_setting: dict[str, list[dict]] = {}
    families_by_setting: dict[str, list[dict]] = {}
    full_by_scale: dict[int, list[dict]] = {}
    for setting, values in settings.items():
        combined = []
        for scale in SCALES:
            windows = sliding_windows(values, scale, setting)
            ordered = sorted(windows, key=lambda row: (-row["mean_score"], row["start"]))
            for rank, window in enumerate(ordered, start=1):
                window["rank_within_scale"] = rank
            combined.extend(windows)
            if setting == "Full":
                full_by_scale[scale] = windows
        windows_by_setting[setting] = combined
        families_by_setting[setting] = cluster_windows(combined)

    top_families = families_by_setting["Full"][:5]
    stability_rows = []
    for candidate_rank, family in enumerate(top_families, start=1):
        candidate_id = f"C{candidate_rank}"
        candidate = family["representative"]
        for setting in SETTINGS:
            match = _best_family_match(candidate, families_by_setting[setting])
            representative = match["representative"]
            stability_rows.append(
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
                    "stable_overlap_gt_0_5": match["overlap"] > 0.5,
                }
            )
    stability = pd.DataFrame(stability_rows)

    scale_matches: dict[str, list[int]] = {}
    for candidate_rank, family in enumerate(top_families, start=1):
        candidate_id = f"C{candidate_rank}"
        candidate = family["representative"]
        stable_scales = []
        for scale in SCALES:
            scale_families = cluster_windows(full_by_scale[scale])[:5]
            if _best_family_match(candidate, scale_families)["overlap"] > 0.5:
                stable_scales.append(scale)
        scale_matches[candidate_id] = stable_scales

    candidates = _candidate_table(table, top_families, stability, scale_matches)
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

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scale in SCALES:
        frame = pd.DataFrame(full_by_scale[scale]).sort_values("rank_within_scale")
        frame.to_csv(args.output_dir / f"window_{scale}_scores.csv", index=False)
    candidates.to_csv(args.output_dir / "candidate_regions_preunblind.csv", index=False)
    stability.to_csv(args.output_dir / "candidate_stability_matrix.csv", index=False)
    _plot_candidate_profile(args.output_dir / "candidate_profile.png", table, candidates)
    _plot_robustness(args.output_dir / "robustness_heatmap.png", stability, candidates)
    rewrite_pdb_bfactors(
        args.pdb,
        args.output_dir / "IL24_181_candidate_score.pdb",
        dict(zip(table["position"].astype(int), table["FinalScore"].astype(float))),
    )
    bepipred = pd.read_csv(args.bepipred)
    _write_final_report(args.report, table, bepipred, candidates, stability)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--feature-table",
        type=Path,
        default=Path("results/07_integration/residue_feature_table.csv"),
    )
    parser.add_argument("--reference", type=Path, default=Path("inputs/reference_181.fasta"))
    parser.add_argument(
        "--alignment",
        type=Path,
        default=Path("results/06_conservation/mafft_alignment.fasta"),
    )
    parser.add_argument("--pdb", type=Path, default=Path("inputs/AF-Q925S4-F1-model_v6.pdb"))
    parser.add_argument(
        "--bepipred",
        type=Path,
        default=Path("results/03_bepipred/bepipred_residue_scores.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/08_candidates"))
    parser.add_argument("--report", type=Path, default=Path("reports/phase1_preunblind_report.md"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
