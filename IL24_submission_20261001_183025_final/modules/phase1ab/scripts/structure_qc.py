#!/usr/bin/env python3
"""Sequence-gated QC for the frozen 181-residue IL-24 AlphaFold PDB."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, median

from Bio import __version__ as BIOPYTHON_VERSION
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


EXPECTED_LENGTH = 181
BACKBONE_ATOMS = ("N", "CA", "C", "O")


class StructureValidationError(RuntimeError):
    """Raised when the PDB cannot be mapped exactly to the frozen reference."""

    def __init__(self, message: str, qc: dict):
        super().__init__(message)
        self.qc = qc


def load_reference(path: str | Path) -> str:
    """Read a single-record FASTA and enforce the frozen project length."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    sequence = sequence.upper()
    if len(sequence) != EXPECTED_LENGTH:
        raise ValueError(
            f"Frozen reference must contain {EXPECTED_LENGTH} residues; found {len(sequence)}"
        )
    return sequence


def _plddt_category(value: float) -> str:
    if value > 90.0:
        return "very_high"
    if value >= 70.0:
        return "high"
    if value >= 50.0:
        return "low"
    return "very_low"


def _discontinuities(residues) -> tuple[list[dict], list[dict]]:
    residue_id_gaps: list[dict] = []
    peptide_bond_gaps: list[dict] = []
    for previous, current in zip(residues, residues[1:]):
        previous_id = previous.id
        current_id = current.id
        sequential = (
            previous_id[0] == " "
            and current_id[0] == " "
            and previous_id[2].strip() == ""
            and current_id[2].strip() == ""
            and current_id[1] == previous_id[1] + 1
        )
        if not sequential:
            residue_id_gaps.append(
                {
                    "after_resseq": previous_id[1],
                    "before_resseq": current_id[1],
                    "after_icode": previous_id[2].strip(),
                    "before_icode": current_id[2].strip(),
                }
            )
        if "C" in previous and "N" in current:
            distance = float(current["N"] - previous["C"])
            if distance > 2.0:
                peptide_bond_gaps.append(
                    {
                        "after_resseq": previous_id[1],
                        "before_resseq": current_id[1],
                        "c_to_n_distance_angstrom": round(distance, 3),
                    }
                )
    return residue_id_gaps, peptide_bond_gaps


def analyze_structure(pdb_path: str | Path, reference: str) -> tuple[list[dict], dict]:
    """Parse the AlphaFold PDB and return project-numbered pLDDT rows plus QC."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("IL24", str(pdb_path))
    models = list(structure.get_models())
    model = models[0]
    chains = list(model.get_chains())
    chain_ids = [chain.id for chain in chains]
    analyzed_chain = chains[0] if len(chains) == 1 else None

    qc: dict = {
        "input_pdb": str(Path(pdb_path).resolve()),
        "reference_length": len(reference),
        "expected_reference_length": EXPECTED_LENGTH,
        "model_count": len(models),
        "chain_ids": chain_ids,
        "analyzed_chain": analyzed_chain.id if analyzed_chain else None,
        "biopython_version": BIOPYTHON_VERSION,
    }
    if len(models) != 1 or analyzed_chain is None or chain_ids != ["A"]:
        qc.update({"analysis_allowed": False, "sequence_exact_match": False})
        raise StructureValidationError(
            "Expected one model containing only analysis chain A", qc
        )

    residues = [residue for residue in analyzed_chain if is_aa(residue, standard=True)]
    pdb_sequence = "".join(seq1(residue.resname, undef_code="X") for residue in residues)
    sequence_exact_match = pdb_sequence == reference
    id_gaps, bond_gaps = _discontinuities(residues)
    missing_backbone: list[dict] = []
    abnormal_residue_ids: list[dict] = []
    rows: list[dict] = []

    for position, residue in enumerate(residues, start=1):
        missing = [atom for atom in BACKBONE_ATOMS if atom not in residue]
        if missing:
            missing_backbone.append(
                {
                    "position": position,
                    "pdb_resseq": residue.id[1],
                    "missing_atoms": missing,
                }
            )
        if residue.id[0] != " " or residue.id[2].strip() or residue.id[1] <= 0:
            abnormal_residue_ids.append(
                {
                    "position": position,
                    "hetflag": residue.id[0].strip(),
                    "pdb_resseq": residue.id[1],
                    "icode": residue.id[2].strip(),
                }
            )
        b_factors = [float(atom.bfactor) for atom in residue.get_atoms()]
        plddt = mean(b_factors) if b_factors else math.nan
        rows.append(
            {
                "position": position,
                "aa": seq1(residue.resname, undef_code="X"),
                "pdb_chain": analyzed_chain.id,
                "pdb_resseq": residue.id[1],
                "pdb_icode": residue.id[2].strip(),
                "plddt": round(plddt, 4),
                "plddt_category": _plddt_category(plddt),
                "uncertainty_flag": plddt < 70.0,
                "backbone_complete": not missing,
            }
        )

    plddt_values = [row["plddt"] for row in rows]
    category_counts = Counter(row["plddt_category"] for row in rows)
    qc.update(
        {
            "residue_count": len(residues),
            "pdb_sequence": pdb_sequence,
            "sequence_exact_match": sequence_exact_match,
            "residue_id_continuous": not id_gaps,
            "residue_id_discontinuities": id_gaps,
            "peptide_bond_discontinuities": bond_gaps,
            "abnormal_residue_ids": abnormal_residue_ids,
            "missing_backbone_atoms": missing_backbone,
            "plddt_summary": {
                "minimum": min(plddt_values) if plddt_values else None,
                "maximum": max(plddt_values) if plddt_values else None,
                "mean": round(mean(plddt_values), 4) if plddt_values else None,
                "median": round(median(plddt_values), 4) if plddt_values else None,
                "category_counts": dict(sorted(category_counts.items())),
                "below_70_count": sum(value < 70.0 for value in plddt_values),
            },
        }
    )
    qc["analysis_allowed"] = bool(
        len(residues) == EXPECTED_LENGTH
        and sequence_exact_match
        and not id_gaps
        and not missing_backbone
        and all(math.isfinite(value) for value in plddt_values)
    )
    if not sequence_exact_match or len(residues) != EXPECTED_LENGTH:
        qc["analysis_allowed"] = False
        raise StructureValidationError(
            "PDB sequence does not exactly match the frozen 181-residue reference", qc
        )
    return rows, qc


def write_structure_outputs(rows: list[dict], qc: dict, output_dir: str | Path) -> dict[str, Path]:
    """Write the Step 2 CSV, JSON, and pLDDT profile."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "residue_plddt.csv"
    json_path = output_dir / "structure_qc.json"
    plot_path = output_dir / "plddt_profile.png"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(
        json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    positions = [row["position"] for row in rows]
    values = [row["plddt"] for row in rows]
    fig, axis = plt.subplots(figsize=(10, 4.2))
    axis.plot(positions, values, color="#275D8C", linewidth=1.5)
    axis.axhline(90, color="#2E8B57", linestyle="--", linewidth=0.8)
    axis.axhline(70, color="#D28E00", linestyle="--", linewidth=0.8)
    axis.axhline(50, color="#B23A48", linestyle="--", linewidth=0.8)
    axis.fill_between(positions, 0, 70, color="#D95F59", alpha=0.08)
    axis.set(xlabel="Project residue position", ylabel="AlphaFold pLDDT", xlim=(1, 181), ylim=(0, 100))
    axis.set_title("IL-24 AlphaFold per-residue confidence")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    return {
        "residue_plddt": csv_path,
        "structure_qc": json_path,
        "plddt_profile": plot_path,
    }


def write_structure_report(qc: dict, report_path: str | Path) -> Path:
    """Write the Step 2 method and result report from computed QC values."""
    summary = qc["plddt_summary"]
    categories = summary["category_counts"]
    text = f"""# Step 2：3D Structure QC

## 1. 本步做了什么

对 AlphaFold PDB 做了序列硬门控、链与残基完整性检查，并从每个 residue 的 B-factor 字段提取 pLDDT。所有结果使用 project 1–181 编号。

## 2. 使用的算法

Biopython PDBParser {qc['biopython_version']}；AlphaFold per-residue pLDDT 读取。

## 3. 算法属于什么

传统结构生物信息学解析与质量控制，不是新的 AI 推理或训练。

## 4. 算法原理

PDBParser 把 PDB 原子记录组织为 model、chain、residue 和 atom 层级。AlphaFold 将局部结构置信度写入 PDB B-factor 字段；本分析取同一 residue 全部原子的均值。pLDDT 仅表示局部坐标可信度，不表示功能重要性。

## 5. 为什么本项目需要它

后续 RSA、DSSP 和结构表位分数必须严格映射到同一条 181-aa 序列；结构低置信位置也需要被显式标记以限制解释。

## 6. 输入

`inputs/AF-Q925S4-F1-model_v6.pdb`；`inputs/reference_181.fasta`。

## 7. 输出

`results/02_structure_qc/residue_plddt.csv`、`structure_qc.json`、`plddt_profile.png`。

## 8. 关键数值结果

- model 数：{qc['model_count']}；chain：{', '.join(qc['chain_ids'])}；残基数：{qc['residue_count']}。
- 与冻结 reference 完全一致：{qc['sequence_exact_match']}；结构分析允许继续：{qc['analysis_allowed']}。
- 缺失主链原子 residue 数：{len(qc['missing_backbone_atoms'])}；residue ID 断点：{len(qc['residue_id_discontinuities'])}；肽键几何断点：{len(qc['peptide_bond_discontinuities'])}。
- pLDDT minimum/mean/median/maximum：{summary['minimum']:.2f} / {summary['mean']:.2f} / {summary['median']:.2f} / {summary['maximum']:.2f}。
- very high (>90)：{categories.get('very_high', 0)}；high (70–90)：{categories.get('high', 0)}；low (50–<70)：{categories.get('low', 0)}；very low (<50)：{categories.get('very_low', 0)}。

## 9. 当前结果怎样解释

pLDDT <70 的位置仅标记结构不确定性。高 pLDDT 不能证明 residue 具有抗原性或功能重要性，pLDDT 不进入最终生物学评分。

## 10. 是否存在问题

序列硬门控状态见 `structure_qc.json`。任何低 pLDDT 位置的结构来源指标应降级解释。

## 11. 下一步

只有序列门控通过后，才使用 FreeSASA 计算 SASA/RSA，并用 DSSP 添加二级结构注释。
"""
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(text, encoding="utf-8")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    reference = load_reference(args.reference)
    try:
        rows, qc = analyze_structure(args.pdb, reference)
    except StructureValidationError as error:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "structure_qc.json").write_text(
            json.dumps(error.qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        raise SystemExit(str(error)) from error
    write_structure_outputs(rows, qc, args.output_dir)
    write_structure_report(qc, args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
