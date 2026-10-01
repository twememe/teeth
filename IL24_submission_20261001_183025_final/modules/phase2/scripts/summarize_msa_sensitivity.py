#!/usr/bin/env python3
"""Compare matched 50-step, three-seed empty-MSA and ColabFold-MSA runs."""

from __future__ import annotations

import csv
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "13_interface_analysis"
OUT = ROOT / "results" / "14_quality"


def load(condition: str) -> list[dict[str, str]]:
    path = BASE / condition / "complex_model_summary.csv"
    return list(csv.DictReader(path.open(encoding="utf-8")))


def epitope(row: dict[str, str]) -> set[int]:
    return {int(value) for value in row["epitope_project_residues"].split(";") if value}


def pairwise_jaccard(rows: list[dict[str, str]]) -> float:
    scores = []
    for left, right in combinations(rows, 2):
        a, b = epitope(left), epitope(right)
        scores.append(len(a & b) / len(a | b) if a | b else 1.0)
    return float(np.mean(scores))


def main() -> None:
    records = {"empty_msa": load("pilot"), "colabfold_msa": load("msa_sensitivity")}
    rows = []
    metrics = ("iptm", "complex_plddt", "hotspot_coverage", "epitope_residues_4p5")
    for condition, data in records.items():
        for route in ("native_blind", "immunogen_blind"):
            selected = [row for row in data if row["route"] == route]
            result: dict[str, object] = {
                "route": route,
                "condition": condition,
                "n_seeds": len(selected),
                "pairwise_epitope_jaccard": round(pairwise_jaccard(selected), 6),
            }
            counts = Counter(residue for row in selected for residue in epitope(row))
            consensus = {residue for residue, count in counts.items() if count >= len(selected) / 2}
            result["epitope_consensus_ge50pct"] = ";".join(map(str, sorted(consensus)))
            for metric in metrics:
                values = np.asarray([float(row[metric]) for row in selected])
                result[f"{metric}_mean"] = round(float(values.mean()), 6)
                result[f"{metric}_sd"] = round(float(values.std(ddof=1)), 6)
            rows.append(result)

    by_key = {(row["route"], row["condition"]): row for row in rows}
    for condition in records:
        native = {int(x) for x in str(by_key[("native_blind", condition)]["epitope_consensus_ge50pct"]).split(";") if x}
        immunogen = {int(x) for x in str(by_key[("immunogen_blind", condition)]["epitope_consensus_ge50pct"]).split(";") if x}
        shared_region = set(range(27, 161))
        native &= shared_region
        immunogen &= shared_region
        cross_jaccard = len(native & immunogen) / len(native | immunogen) if native | immunogen else 1.0
        for route in ("native_blind", "immunogen_blind"):
            by_key[(route, condition)]["native_immunogen_consensus_jaccard"] = round(cross_jaccard, 6)
    for route in ("native_blind", "immunogen_blind"):
        empty = by_key[(route, "empty_msa")]
        msa = by_key[(route, "colabfold_msa")]
        for metric in ("pairwise_epitope_jaccard",) + tuple(f"{name}_mean" for name in metrics):
            msa[f"delta_vs_empty_{metric}"] = round(float(msa[metric]) - float(empty[metric]), 6)

    OUT.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with (OUT / "msa_sensitivity_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} matched-condition summaries")


if __name__ == "__main__":
    main()
