#!/usr/bin/env python3
"""Generate the standardized Phase 2 candidate list from completed analyses."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RANK = ROOT / "results/14_quality/model_ranking_formal_msa.csv"
if not RANK.exists():
    RANK = ROOT / "results/model_ranking_formal_msa.csv"
OUT = ROOT / "results.csv"

def main():
    rows = list(csv.DictReader(RANK.open(encoding="utf-8")))
    required = {"route", "seed", "confidence_score", "iptm", "complex_plddt", "prodigy_dG_kcal_mol"}
    if not rows:
        raise ValueError(f"排序文件为空: {RANK}")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"缺少字段: {sorted(missing)}")
    for r in rows:
        try:
            float(r["confidence_score"]); float(r["iptm"]); float(r["complex_plddt"]); float(r["prodigy_dG_kcal_mol"])
        except (TypeError, ValueError) as e:
            raise ValueError(f"模型数值字段无法解析: {r}") from e
    rows.sort(key=lambda r: (float(r["confidence_score"]), float(r["iptm"])), reverse=True)
    out = []
    for i, r in enumerate(rows[:10], 1):
        route, seed = r["route"], r["seed"]
        full_structure = ROOT / "results/12_complex_prediction/formal_msa" / f"{route}_seed{seed}" / f"boltz_results_{route}_msa/predictions/{route}_msa/{route}_msa_model_0.pdb"
        if full_structure.exists():
            structure = str(full_structure.relative_to(ROOT)).replace("\\", "/")
        else:
            # Compact release ships representative structures rather than all raw predictions.
            structure = f"results/{route}_seed{seed}.pdb"
        if not (ROOT / structure).is_file():
            raise FileNotFoundError(f"结构文件不存在: {structure}")
        out.append({
            "candidate_id": f"IA6-13-8-{i:02d}",
            "track": "赛道一：AI大分子与多肽药物设计",
            "task": "IL-24 antibody structure prediction",
            "route": route,
            "seed": seed,
            "msa_condition": "formal_msa",
            "sequence_or_structure": structure,
            "source_pdb": structure,
            "structure_status": "unrelaxed",
            "relaxation_status": "not_relaxed",
            "clash_qc": "see_relaxed_geometry_qc; raw model may contain severe clashes",
            "model_version": "Boltz-1 2.2.1; formal MSA; seed " + seed,
            "ipTM": r["iptm"],
            "complex_pLDDT": r["complex_plddt"],
            "PRODIGY_dG_kcal_mol": r["prodigy_dG_kcal_mol"],
            "hotspot_coverage": r["hotspot_coverage"],
            "ranking_basis": f"confidence_rank={r['confidence_rank']}; consensus_rank={r['consensus_rank']}",
            "notes": "Computational ranking only; PRODIGY is relative and not experimental KD."
        })
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    print(f"Wrote {OUT} with {len(out)} candidates")

if __name__ == "__main__":
    main()
