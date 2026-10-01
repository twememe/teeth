#!/usr/bin/env python
"""Integrate the three available IL-24 Phase 1 evidence groups."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


FEATURE_COLUMNS = [
    "position",
    "AA",
    "pLDDT",
    "BepiPred_raw",
    "BepiPred_percentile",
    "SASA",
    "RSA_raw",
    "RSA_percentile",
    "Conservation_raw",
    "Conservation_percentile",
    "SecondaryStructure",
    "FinalScore",
    "StructureConfidenceFlag",
]


def read_fasta(path: Path) -> str:
    sequence = "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )
    if not sequence:
        raise ValueError(f"No sequence in {path}")
    return sequence


def percentile_rank(values: list[float]) -> list[float]:
    """Average-tie ranks scaled to exact 0 and 1 endpoints."""
    numbers = pd.Series(values, dtype=float)
    if not numbers.map(math.isfinite).all():
        raise ValueError("Percentile input contains NaN or infinity")
    if len(numbers) == 1:
        return [0.5]
    ranks = numbers.rank(method="average")
    return ((ranks - 1.0) / (len(numbers) - 1.0)).tolist()


def _validate_mapping(reference: str, frame: pd.DataFrame, aa_column: str, label: str) -> None:
    if len(frame) != len(reference):
        raise ValueError(f"{label} row count {len(frame)} != reference length {len(reference)}")
    if frame["position"].astype(int).tolist() != list(range(1, len(reference) + 1)):
        raise ValueError(f"{label} positions are not continuous project numbering")
    observed = "".join(frame[aa_column].astype(str).str.upper())
    if observed != reference:
        mismatch = next(index for index, pair in enumerate(zip(observed, reference), 1) if pair[0] != pair[1])
        raise ValueError(f"AA mismatch in {label} at project position {mismatch}")


def _confidence_flag(plddt: float) -> str:
    if plddt >= 90.0:
        return "VERY_HIGH"
    if plddt >= 70.0:
        return "HIGH"
    if plddt >= 50.0:
        return "LOW"
    return "VERY_LOW"


def assemble_feature_table(
    reference: str,
    bepipred: pd.DataFrame,
    plddt: pd.DataFrame,
    sasa: pd.DataFrame,
    conservation: pd.DataFrame,
    secondary: pd.DataFrame,
) -> pd.DataFrame:
    _validate_mapping(reference, bepipred, "AA", "BepiPred")
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
            "SASA": pd.to_numeric(sasa["absolute_sasa_angstrom2"], errors="raise"),
            "RSA_raw": pd.to_numeric(sasa["rsa"], errors="raise"),
            "Conservation_raw": pd.to_numeric(
                conservation["conservation_score"], errors="raise"
            ),
            "SecondaryStructure": secondary["secondary_structure_class"].astype(str),
        }
    )
    for column in ["pLDDT", "BepiPred_raw", "SASA", "RSA_raw", "Conservation_raw"]:
        if not table[column].map(math.isfinite).all():
            raise ValueError(f"Non-finite values in {column}")
    table["BepiPred_percentile"] = percentile_rank(table["BepiPred_raw"].tolist())
    table["RSA_percentile"] = percentile_rank(table["RSA_raw"].tolist())
    table["Conservation_percentile"] = percentile_rank(table["Conservation_raw"].tolist())
    table["FinalScore"] = table[
        ["BepiPred_percentile", "RSA_percentile", "Conservation_percentile"]
    ].mean(axis=1)
    table["StructureConfidenceFlag"] = table["pLDDT"].map(_confidence_flag)
    return table[FEATURE_COLUMNS]


def _plot_correlation(path: Path, table: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    columns = [
        "BepiPred_percentile",
        "RSA_percentile",
        "Conservation_percentile",
        "FinalScore",
    ]
    labels = ["BepiPred", "RSA", "Conservation", "Final"]
    correlation = table[columns].corr().to_numpy()
    figure, axis = plt.subplots(figsize=(6.2, 5.4))
    image = axis.imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
    axis.set_xticks(range(4), labels=labels, rotation=35, ha="right")
    axis.set_yticks(range(4), labels=labels)
    for row in range(4):
        for column in range(4):
            axis.text(column, row, f"{correlation[row, column]:.2f}", ha="center", va="center")
    axis.set_title("IL-24 available-evidence correlation")
    figure.colorbar(image, ax=axis, shrink=0.82, label="Pearson r")
    figure.tight_layout()
    figure.savefig(path, dpi=220)
    plt.close(figure)


def _plot_tracks(path: Path, table: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = table["position"]
    figure, axes = plt.subplots(4, 1, figsize=(12, 8.5), sharex=True)
    tracks = [
        ("BepiPred_percentile", "#2166ac", "BepiPred percentile"),
        ("RSA_percentile", "#1b9e77", "RSA percentile"),
        ("Conservation_percentile", "#d95f02", "Conservation percentile"),
        ("FinalScore", "#762a83", "Available-evidence FinalScore"),
    ]
    for axis, (column, color, label) in zip(axes, tracks):
        axis.plot(x, table[column], color=color, linewidth=1.25)
        axis.set_ylabel(label)
        axis.set_ylim(-0.03, 1.03)
        axis.grid(alpha=0.18)
    axes[-1].set_xlabel("Project position (1–181)")
    axes[-1].set_xlim(1, len(table))
    figure.suptitle("IL-24 Phase 1 feature tracks (DiscoTope unavailable)")
    figure.tight_layout()
    figure.savefig(path, dpi=220)
    plt.close(figure)


def _write_report(path: Path, table: pd.DataFrame) -> None:
    correlations = table[
        ["BepiPred_percentile", "RSA_percentile", "Conservation_percentile"]
    ].corr()
    path.write_text(
        f"""# Step 7 — Multi-Evidence Integration

## Status

`discotope_available = false`. The preregistered DiscoTope-3.0 structural-epitope branch is missing because both allowed official installation attempts ended in external dependency retrieval failure. No substitute model or fabricated 3D epitope score was introduced.

## Available-evidence consensus

Each available continuous feature was independently transformed with average-tie percentile ranks scaled to [0, 1]: BepiPred-3.0 probability (B), FreeSASA RSA (R), and normalized Shannon conservation (C).

`FinalScore = (B + R + C) / 3`

pLDDT is retained only as a structure-confidence flag. DSSP is retained only as annotation; neither contributes to FinalScore.

## Validation and summary

- Rows: {len(table)}; project positions: 1–{int(table['position'].max())}.
- FinalScore range: {table['FinalScore'].min():.6f}–{table['FinalScore'].max():.6f}; mean {table['FinalScore'].mean():.6f}.
- BepiPred–RSA Pearson r: {correlations.loc['BepiPred_percentile', 'RSA_percentile']:.4f}.
- BepiPred–Conservation Pearson r: {correlations.loc['BepiPred_percentile', 'Conservation_percentile']:.4f}.
- RSA–Conservation Pearson r: {correlations.loc['RSA_percentile', 'Conservation_percentile']:.4f}.

These are multi-evidence computational prioritization scores, not experimentally validated epitopes, receptor interfaces, or confirmed neutralizing regions.
""",
        encoding="utf-8",
    )


def run(args: argparse.Namespace) -> None:
    reference = read_fasta(args.reference)
    if len(reference) != 181:
        raise ValueError(f"Frozen reference must contain exactly 181 aa; found {len(reference)}")
    table = assemble_feature_table(
        reference,
        pd.read_csv(args.bepipred),
        pd.read_csv(args.plddt),
        pd.read_csv(args.sasa),
        pd.read_csv(args.conservation),
        pd.read_csv(args.secondary),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "residue_feature_table.csv", index=False)
    consensus = table[
        [
            "position",
            "AA",
            "BepiPred_percentile",
            "RSA_percentile",
            "Conservation_percentile",
            "FinalScore",
            "StructureConfidenceFlag",
        ]
    ].copy()
    consensus["discotope_available"] = False
    consensus.to_csv(args.output_dir / "residue_consensus_scores.csv", index=False)
    _plot_correlation(args.output_dir / "feature_correlation.png", table)
    _plot_tracks(args.output_dir / "all_feature_tracks.png", table)
    _write_report(args.report, table)
    (args.output_dir / "provenance.json").write_text(
        json.dumps(
            {
                "discotope_available": False,
                "discotope_status": "MISSING_EXTERNAL_DEPENDENCY_FAILURE",
                "primary_consensus": "available_evidence_equal_weight",
                "evidence_groups": ["BepiPred", "RSA", "Conservation"],
                "formula": "(BepiPred_percentile + RSA_percentile + Conservation_percentile) / 3",
                "pLDDT_role": "structure_confidence_flag_only",
                "DSSP_role": "annotation_only",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("inputs/reference_181.fasta"))
    parser.add_argument(
        "--bepipred",
        type=Path,
        default=Path("results/03_bepipred/bepipred_residue_scores.csv"),
    )
    parser.add_argument("--plddt", type=Path, default=Path("results/02_structure_qc/residue_plddt.csv"))
    parser.add_argument("--sasa", type=Path, default=Path("results/04_surface/residue_sasa.csv"))
    parser.add_argument(
        "--secondary",
        type=Path,
        default=Path("results/04_surface/residue_secondary_structure.csv"),
    )
    parser.add_argument(
        "--conservation",
        type=Path,
        default=Path("results/06_conservation/residue_conservation.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/07_integration"))
    parser.add_argument("--report", type=Path, default=Path("reports/step07_integration.md"))
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
