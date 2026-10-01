#!/usr/bin/env python3
"""Integrate restored DiscoTope evidence without mutating frozen Phase 1A."""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import pandas as pd

try:
    from scripts.integration import (
        _confidence_flag,
        _validate_mapping,
        percentile_rank,
        read_fasta,
    )
except ModuleNotFoundError:
    from integration import (
        _confidence_flag,
        _validate_mapping,
        percentile_rank,
        read_fasta,
    )


FEATURE_COLUMNS = [
    "position",
    "AA",
    "pLDDT",
    "BepiPred_raw",
    "BepiPred_percentile",
    "DiscoTope_raw",
    "DiscoTope_percentile",
    "EpitopeConsensus",
    "SASA",
    "RSA_raw",
    "RSA_percentile",
    "Conservation_raw",
    "Conservation_percentile",
    "SecondaryStructure",
    "FinalScore_Phase1A",
    "FinalScore_Phase1B",
    "FinalScore_4way_sensitivity",
    "StructureConfidenceFlag",
]


def assemble_phase1b_feature_table(
    reference: str,
    bepipred: pd.DataFrame,
    discotope: pd.DataFrame,
    plddt: pd.DataFrame,
    sasa: pd.DataFrame,
    conservation: pd.DataFrame,
    secondary: pd.DataFrame,
) -> pd.DataFrame:
    _validate_mapping(reference, bepipred, "AA", "BepiPred")
    _validate_mapping(reference, discotope, "AA", "DiscoTope")
    _validate_mapping(reference, plddt, "aa", "pLDDT")
    _validate_mapping(reference, sasa, "aa", "SASA/RSA")
    _validate_mapping(reference, conservation, "AA", "conservation")
    _validate_mapping(reference, secondary, "aa", "secondary structure")

    table = pd.DataFrame(
        {
            "position": range(1, len(reference) + 1),
            "AA": list(reference),
            "pLDDT": pd.to_numeric(plddt["plddt"], errors="raise"),
            "BepiPred_raw": pd.to_numeric(
                bepipred["bepipred_positive_probability"], errors="raise"
            ),
            "DiscoTope_raw": pd.to_numeric(
                discotope["DiscoTope_raw"], errors="raise"
            ),
            "SASA": pd.to_numeric(
                sasa["absolute_sasa_angstrom2"], errors="raise"
            ),
            "RSA_raw": pd.to_numeric(sasa["rsa"], errors="raise"),
            "Conservation_raw": pd.to_numeric(
                conservation["conservation_score"], errors="raise"
            ),
            "SecondaryStructure": secondary[
                "secondary_structure_class"
            ].astype(str),
        }
    )
    numeric = [
        "pLDDT",
        "BepiPred_raw",
        "DiscoTope_raw",
        "SASA",
        "RSA_raw",
        "Conservation_raw",
    ]
    for column in numeric:
        if not table[column].map(math.isfinite).all():
            raise ValueError(f"Non-finite values in {column}")

    table["BepiPred_percentile"] = percentile_rank(
        table["BepiPred_raw"].tolist()
    )
    table["DiscoTope_percentile"] = percentile_rank(
        table["DiscoTope_raw"].tolist()
    )
    table["RSA_percentile"] = percentile_rank(table["RSA_raw"].tolist())
    table["Conservation_percentile"] = percentile_rank(
        table["Conservation_raw"].tolist()
    )
    table["EpitopeConsensus"] = table[
        ["BepiPred_percentile", "DiscoTope_percentile"]
    ].mean(axis=1)
    table["FinalScore_Phase1A"] = table[
        ["BepiPred_percentile", "RSA_percentile", "Conservation_percentile"]
    ].mean(axis=1)
    table["FinalScore_Phase1B"] = table[
        ["EpitopeConsensus", "RSA_percentile", "Conservation_percentile"]
    ].mean(axis=1)
    table["FinalScore_4way_sensitivity"] = table[
        [
            "BepiPred_percentile",
            "DiscoTope_percentile",
            "RSA_percentile",
            "Conservation_percentile",
        ]
    ].mean(axis=1)
    table["StructureConfidenceFlag"] = table["pLDDT"].map(_confidence_flag)
    return table[FEATURE_COLUMNS]


def correlation_summary(table: pd.DataFrame) -> dict[str, float]:
    return {
        "pearson_raw": float(
            table["BepiPred_raw"].corr(table["DiscoTope_raw"], method="pearson")
        ),
        "spearman_raw": float(
            table["BepiPred_raw"].corr(table["DiscoTope_raw"], method="spearman")
        ),
        "pearson_percentile": float(
            table["BepiPred_percentile"].corr(
                table["DiscoTope_percentile"], method="pearson"
            )
        ),
        "spearman_percentile": float(
            table["BepiPred_percentile"].corr(
                table["DiscoTope_percentile"], method="spearman"
            )
        ),
    }


def _plot_scatter(path: Path, table: pd.DataFrame, correlations: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    x = table["BepiPred_raw"].to_numpy()
    y = table["DiscoTope_raw"].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    fit_x = np.linspace(x.min(), x.max(), 100)
    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    points = axis.scatter(
        x,
        y,
        c=table["pLDDT"],
        cmap="viridis",
        s=32,
        alpha=0.78,
        edgecolor="none",
    )
    axis.plot(fit_x, slope * fit_x + intercept, color="#C62828", linewidth=1.4)
    axis.text(
        0.03,
        0.97,
        f"Pearson r = {correlations['pearson_raw']:.3f}\n"
        f"Spearman rho = {correlations['spearman_raw']:.3f}",
        transform=axis.transAxes,
        va="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )
    axis.set_xlabel("BepiPred-3.0 positive probability")
    axis.set_ylabel("DiscoTope-3.0 official raw score")
    axis.set_title("Sequence vs conformational epitope evidence")
    axis.grid(alpha=0.18)
    figure.colorbar(points, ax=axis, label="pLDDT (annotation only)")
    figure.tight_layout()
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def _plot_tracks(path: Path, table: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = table["position"]
    tracks = [
        ("BepiPred_percentile", "#1565C0", "BepiPred percentile"),
        ("DiscoTope_percentile", "#6A1B9A", "DiscoTope percentile"),
        ("RSA_percentile", "#00897B", "RSA percentile"),
        ("Conservation_percentile", "#EF6C00", "Conservation percentile"),
        ("FinalScore_Phase1B", "#C62828", "Primary FinalScore Phase 1B"),
    ]
    figure, axes = plt.subplots(5, 1, figsize=(13, 10), sharex=True)
    for axis, (column, color, label) in zip(axes, tracks):
        axis.plot(x, table[column], color=color, linewidth=1.35)
        axis.fill_between(x, 0, table[column], color=color, alpha=0.12)
        axis.set_ylabel(label)
        axis.set_ylim(-0.03, 1.03)
        axis.grid(axis="y", alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Project position (1-181)")
    axes[-1].set_xlim(1, 181)
    figure.suptitle("IL-24 Phase 1B complete evidence tracks")
    figure.tight_layout()
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def _plot_phase_comparison(path: Path, table: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = table["position"]
    delta = table["FinalScore_Phase1B"] - table["FinalScore_Phase1A"]
    figure, axes = plt.subplots(
        2, 1, figsize=(13, 7.2), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )
    axes[0].plot(x, table["FinalScore_Phase1A"], color="#607D8B", label="Phase 1A")
    axes[0].plot(x, table["FinalScore_Phase1B"], color="#C62828", label="Phase 1B")
    axes[0].set_ylabel("Residue consensus score")
    axes[0].legend(frameon=False, ncol=2)
    axes[0].grid(axis="y", alpha=0.2)
    colors = ["#2E7D32" if value >= 0 else "#C62828" for value in delta]
    axes[1].bar(x, delta, color=colors, width=1.0)
    axes[1].axhline(0, color="#263238", linewidth=0.8)
    axes[1].set_ylabel("Phase 1B - 1A")
    axes[1].set_xlabel("Project position (1-181)")
    axes[1].set_xlim(1, 181)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("Residue-level effect of restored 3D epitope evidence")
    figure.tight_layout()
    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def run(args: argparse.Namespace) -> None:
    reference = read_fasta(args.reference)
    if len(reference) != 181:
        raise ValueError(f"Frozen reference must be 181 aa; found {len(reference)}")
    table = assemble_phase1b_feature_table(
        reference,
        pd.read_csv(args.bepipred),
        pd.read_csv(args.discotope),
        pd.read_csv(args.plddt),
        pd.read_csv(args.sasa),
        pd.read_csv(args.conservation),
        pd.read_csv(args.secondary),
    )
    frozen = pd.read_csv(args.phase1a_table)
    _validate_mapping(reference, frozen, "AA", "frozen Phase 1A")
    reproduction_error = float(
        (table["FinalScore_Phase1A"] - frozen["FinalScore"]).abs().max()
    )
    if reproduction_error >= 1e-12:
        raise RuntimeError(
            f"BepiPred-only score failed Phase 1A reproduction: {reproduction_error}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "residue_feature_table_3d.csv"
    table.to_csv(table_path, index=False)
    table[
        [
            "position",
            "AA",
            "BepiPred_percentile",
            "DiscoTope_percentile",
            "EpitopeConsensus",
            "RSA_percentile",
            "Conservation_percentile",
            "FinalScore_Phase1A",
            "FinalScore_Phase1B",
            "FinalScore_4way_sensitivity",
            "StructureConfidenceFlag",
        ]
    ].to_csv(args.output_dir / "residue_consensus_scores_3d.csv", index=False)

    delta = table[
        [
            "position",
            "AA",
            "BepiPred_percentile",
            "DiscoTope_percentile",
            "EpitopeConsensus",
            "FinalScore_Phase1A",
            "FinalScore_Phase1B",
        ]
    ].copy()
    delta["Delta_FinalScore_B_minus_A"] = (
        delta["FinalScore_Phase1B"] - delta["FinalScore_Phase1A"]
    )
    delta.to_csv(
        args.output_dir / "phase1a_vs_phase1b_residue_delta.csv", index=False
    )

    correlations = correlation_summary(table)
    (args.output_dir / "bepipred_discotope_correlations.json").write_text(
        json.dumps(correlations, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    scatter = args.output_dir / "BepiPred_vs_DiscoTope_scatter.png"
    _plot_scatter(scatter, table, correlations)
    shutil.copyfile(scatter, args.output_dir / "04_bepipred_vs_discotope.png")
    tracks = args.output_dir / "01_all_evidence_tracks_phase1b.png"
    _plot_tracks(tracks, table)
    comparison = args.output_dir / "phase1a_vs_phase1b_profile.png"
    _plot_phase_comparison(comparison, table)
    shutil.copyfile(comparison, args.output_dir / "05_phase1a_vs_phase1b_profile.png")

    delta_sorted_up = delta.sort_values(
        ["Delta_FinalScore_B_minus_A", "position"], ascending=[False, True]
    ).head(10)
    delta_sorted_down = delta.sort_values(
        ["Delta_FinalScore_B_minus_A", "position"], ascending=[True, True]
    ).head(10)
    serious_conflict = int(
        (
            (table["BepiPred_percentile"] >= 0.8)
            & (table["DiscoTope_percentile"] <= 0.2)
        ).sum()
        + (
            (table["BepiPred_percentile"] <= 0.2)
            & (table["DiscoTope_percentile"] >= 0.8)
        ).sum()
    )
    args.report.write_text(
        f"""# Step 7B — Complete integration with restored 3D evidence

## Primary formula

Each raw feature uses the frozen Phase 1A average-tie percentile transform on [0,1]. `DiscoTope_percentile` is computed from the official uncalibrated `DiscoTope-3.0_score` preserved as `DiscoTope_raw`.

`EpitopeConsensus = (BepiPred_percentile + DiscoTope_percentile) / 2`

`FinalScore_Phase1B = (EpitopeConsensus + RSA_percentile + Conservation_percentile) / 3`

The direct four-feature mean is retained only as `FinalScore_4way_sensitivity` and is not used for primary candidate ranking. pLDDT is a confidence flag only; DSSP is annotation only.

## Validation

- Rows: 181; reference/numbering exact.
- BepiPred-only reproduction of frozen Phase 1A: maximum absolute error {reproduction_error:.3e}.
- Primary Phase 1B score range: {table['FinalScore_Phase1B'].min():.6f}–{table['FinalScore_Phase1B'].max():.6f}.

## BepiPred versus DiscoTope

- Raw Pearson r: {correlations['pearson_raw']:.4f}.
- Raw Spearman rho: {correlations['spearman_raw']:.4f}.
- Percentile Pearson r: {correlations['pearson_percentile']:.4f}.
- Percentile Spearman rho: {correlations['spearman_percentile']:.4f}.
- Extreme cross-model conflict residues (one >=0.8 while the other <=0.2): {serious_conflict}/181.

Low or moderate correlation is not a failure: the sequence model and conformational structure model interrogate complementary epitope properties.

## Phase 1A versus Phase 1B residue changes

- Delta range: {delta['Delta_FinalScore_B_minus_A'].min():.6f} to {delta['Delta_FinalScore_B_minus_A'].max():.6f}; mean {delta['Delta_FinalScore_B_minus_A'].mean():.6f}.
- Largest positive deltas: {', '.join(f"{row.AA}{int(row.position)} ({row.Delta_FinalScore_B_minus_A:+.4f})" for row in delta_sorted_up.itertuples())}.
- Largest negative deltas: {', '.join(f"{row.AA}{int(row.position)} ({row.Delta_FinalScore_B_minus_A:+.4f})" for row in delta_sorted_down.itertuples())}.

No historical experimental interval was accessed. These are computational score comparisons only.
""",
        encoding="utf-8",
    )
    (args.output_dir / "provenance.json").write_text(
        json.dumps(
            {
                "phase": "Phase 1B post-freeze user-authorized restoration",
                "historical_region_still_blinded": True,
                "discotope_available": True,
                "percentile_method": "average ties; (rank-1)/(n-1)",
                "primary_formula": "(((B+D)/2)+R+C)/3",
                "sensitivity_4way": "(B+D+R+C)/4; SENSITIVITY_ONLY",
                "discoTope_percentile_source": "official DiscoTope-3.0_score",
                "phase1a_reproduction_max_abs_error": reproduction_error,
                "correlations": correlations,
                "pLDDT_role": "structure confidence flag only",
                "DSSP_role": "annotation only",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("inputs/reference_181.fasta"))
    parser.add_argument("--bepipred", type=Path, default=Path("results/03_bepipred/bepipred_residue_scores.csv"))
    parser.add_argument("--discotope", type=Path, default=Path("results/05b_discotope_restored/discotope_residue_scores.csv"))
    parser.add_argument("--plddt", type=Path, default=Path("results/02_structure_qc/residue_plddt.csv"))
    parser.add_argument("--sasa", type=Path, default=Path("results/04_surface/residue_sasa.csv"))
    parser.add_argument("--conservation", type=Path, default=Path("results/06_conservation/residue_conservation.csv"))
    parser.add_argument("--secondary", type=Path, default=Path("results/04_surface/residue_secondary_structure.csv"))
    parser.add_argument("--phase1a-table", type=Path, default=Path("results/07_integration/residue_feature_table.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/07b_integration_with_3d"))
    parser.add_argument("--report", type=Path, default=Path("reports/step07b_integration_with_3d.md"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
