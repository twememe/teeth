#!/usr/bin/env python
"""Sequence QC for the frozen 181-aa IL-24 project reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
from dataclasses import dataclass
from pathlib import Path

import Bio
from Bio import Align
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class ExtractedChain:
    chain_id: str
    all_chain_ids: list[str]
    sequence: str
    project_positions: list[int]
    pdb_residue_ids: list[str]
    pdb_residue_numbers: list[int]
    pdb_insertion_codes: list[str]
    missing_project_positions: list[int]
    missing_backbone_atoms: dict[int, list[str]]


def read_fasta(path: Path) -> FastaRecord:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header = None
        sequence_parts: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"Expected one FASTA record in {path}")
                header = line[1:]
            elif header is None:
                raise ValueError(f"FASTA sequence precedes header in {path}")
            else:
                sequence_parts.append(line)
    if header is None or not sequence_parts:
        raise ValueError(f"No FASTA record found in {path}")
    sequence = "".join(sequence_parts).replace(" ", "").upper()
    if not sequence.isalpha():
        raise ValueError(f"Non-letter character found in sequence from {path}")
    identifier, _, description = header.partition(" ")
    return FastaRecord(identifier, description, sequence)


def _format_residue_id(number: int, insertion_code: str) -> str:
    return f"{number}{insertion_code}" if insertion_code else str(number)


def extract_chain(path: Path, expected_chain: str = "A") -> ExtractedChain:
    structure = PDBParser(QUIET=True).get_structure("IL24", str(path))
    model = structure[0]
    protein_chains = [
        chain for chain in model if any(is_aa(residue, standard=True) for residue in chain)
    ]
    all_chain_ids = [chain.id for chain in protein_chains]
    if expected_chain not in all_chain_ids:
        raise ValueError(f"Expected protein chain {expected_chain!r}; found {all_chain_ids}")
    chain = model[expected_chain]
    residues = [residue for residue in chain if is_aa(residue, standard=True)]

    sequence = "".join(seq1(residue.resname, undef_code="X") for residue in residues)
    numbers = [int(residue.id[1]) for residue in residues]
    insertion_codes = [str(residue.id[2]).strip() for residue in residues]
    residue_ids = [
        _format_residue_id(number, insertion_code)
        for number, insertion_code in zip(numbers, insertion_codes)
    ]

    missing_project_positions: list[int] = []
    for index, (left, right) in enumerate(zip(numbers, numbers[1:]), start=1):
        if right > left + 1:
            missing_project_positions.extend(range(index + 1, index + (right - left)))

    missing_backbone_atoms: dict[int, list[str]] = {}
    for project_position, residue in enumerate(residues, start=1):
        missing = [atom for atom in ("N", "CA", "C") if atom not in residue]
        if missing:
            missing_backbone_atoms[project_position] = missing

    return ExtractedChain(
        chain_id=expected_chain,
        all_chain_ids=all_chain_ids,
        sequence=sequence,
        project_positions=list(range(1, len(residues) + 1)),
        pdb_residue_ids=residue_ids,
        pdb_residue_numbers=numbers,
        pdb_insertion_codes=insertion_codes,
        missing_project_positions=missing_project_positions,
        missing_backbone_atoms=missing_backbone_atoms,
    )


def align_to_reference(reference: str, subject: str, subject_name: str) -> tuple[list[dict], dict]:
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(reference, subject)[0]
    reference_blocks, subject_blocks = alignment.aligned

    rows: list[dict] = []
    ref_cursor = 0
    subject_cursor = 0
    matches = substitutions = insertions = deletions = 0

    def append_gap_rows(next_ref: int, next_subject: int) -> None:
        nonlocal ref_cursor, subject_cursor, insertions, deletions
        while ref_cursor < next_ref:
            rows.append(
                {
                    "project_position": ref_cursor + 1,
                    "project_aa": reference[ref_cursor],
                    f"{subject_name}_position": "",
                    f"{subject_name}_aa": "",
                    "relation": "deletion",
                    "position_offset": "",
                }
            )
            ref_cursor += 1
            deletions += 1
        while subject_cursor < next_subject:
            rows.append(
                {
                    "project_position": "",
                    "project_aa": "",
                    f"{subject_name}_position": subject_cursor + 1,
                    f"{subject_name}_aa": subject[subject_cursor],
                    "relation": "insertion",
                    "position_offset": "",
                }
            )
            subject_cursor += 1
            insertions += 1

    for (ref_start, ref_end), (subject_start, subject_end) in zip(
        reference_blocks, subject_blocks
    ):
        append_gap_rows(int(ref_start), int(subject_start))
        block_length = int(ref_end - ref_start)
        if block_length != int(subject_end - subject_start):
            raise RuntimeError("Unexpected unequal aligned block lengths")
        for _ in range(block_length):
            relation = "match" if reference[ref_cursor] == subject[subject_cursor] else "substitution"
            matches += relation == "match"
            substitutions += relation == "substitution"
            rows.append(
                {
                    "project_position": ref_cursor + 1,
                    "project_aa": reference[ref_cursor],
                    f"{subject_name}_position": subject_cursor + 1,
                    f"{subject_name}_aa": subject[subject_cursor],
                    "relation": relation,
                    "position_offset": subject_cursor - ref_cursor,
                }
            )
            ref_cursor += 1
            subject_cursor += 1
    append_gap_rows(len(reference), len(subject))

    offsets = sorted(
        {row["position_offset"] for row in rows if row["position_offset"] != ""}
    )
    summary = {
        "subject": subject_name,
        "reference_length": len(reference),
        "subject_length": len(subject),
        "alignment_score": float(alignment.score),
        "match_count": matches,
        "substitution_count": substitutions,
        "insertion_count": insertions,
        "deletion_count": deletions,
        "reference_match_fraction": matches / len(reference),
        "exact_match": reference == subject,
        "numbering_offsets": offsets,
        "substitution_project_positions": [
            row["project_position"] for row in rows if row["relation"] == "substitution"
        ],
        "deleted_project_positions": [
            row["project_position"] for row in rows if row["relation"] == "deletion"
        ],
        "inserted_subject_positions": [
            row[f"{subject_name}_position"] for row in rows if row["relation"] == "insertion"
        ],
    }
    return rows, summary


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_fasta(path: Path, identifier: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{identifier}\n")
        for start in range(0, len(sequence), 60):
            handle.write(sequence[start : start + 60] + "\n")


def _write_report(path: Path, qc: dict, pdb_summary: dict, uniprot_summary: dict) -> None:
    pdb_allowed = qc["pdb_structure_linkage_allowed"]
    substitutions = uniprot_summary["substitution_project_positions"]
    report = f"""# Step 1 — Sequence QC

## Step 1 完成

### 1. 本步做了什么
从 AlphaFold PDB 的 Chain A 提取氨基酸序列，并把 PDB、UniProt Q925S4 与冻结的 181-aa 项目参考逐位对齐。UniProt 只用于映射和来源注释，未替换项目坐标。

### 2. 使用的算法
Biopython PDBParser 与 Biopython PairwiseAligner（Biopython {Bio.__version__}）。

### 3. 算法属于什么
传统生物信息学的结构文件解析和全局双序列比对，不是新的深度学习训练。

### 4. 算法原理
通俗地说，它像把两篇文章逐字排齐，记录相同字符、替换、插入和删除。技术上使用全局 PairwiseAligner，match=2、mismatch=-1、gap open=-5、gap extend=-1，并用 PDBParser 读取标准氨基酸 residue 与 N/CA/C 主链原子。

### 5. 为什么本项目需要它
所有后续“第 N 位”分数都依赖同一坐标；此步骤防止 PDB、项目参考与 UniProt 编号被混用。

### 6. 输入
- `inputs/reference_181.fasta`
- `inputs/AF-Q925S4-F1-model_v6.pdb`
- `inputs/uniprotkb_accession_Q925S4_2026_08_20.fasta.gz`

### 7. 输出
- `results/01_sequence_qc/sequence_mapping.csv`
- `results/01_sequence_qc/pdb_sequence.fasta`
- `results/01_sequence_qc/uniprot_mapping.csv`
- `results/01_sequence_qc/sequence_qc.json`
- `reports/step01_sequence_qc.md`

### 8. 关键数值结果
- Project reference 长度：{qc['reference_length']} aa。
- PDB protein chains：{', '.join(qc['pdb_chain_ids'])}；Chain A 长度：{qc['pdb_sequence_length']} aa。
- Reference↔PDB：matches={pdb_summary['match_count']}，substitutions={pdb_summary['substitution_count']}，insertions={pdb_summary['insertion_count']}，deletions={pdb_summary['deletion_count']}，reference match fraction={pdb_summary['reference_match_fraction']:.6f}，PDB residue-number offset(s)={qc['pdb_residue_number_offsets']}。
- Reference↔UniProt：UniProt 长度={uniprot_summary['subject_length']} aa，matches={uniprot_summary['match_count']}，substitutions={uniprot_summary['substitution_count']}，insertions={uniprot_summary['insertion_count']}，deletions={uniprot_summary['deletion_count']}，offset(s)={uniprot_summary['numbering_offsets']}，substitution project position(s)={substitutions}。
- PDB 结构联动允许：{pdb_allowed}。

### 9. 当前结果怎样解释
本结果只证明序列对象和编号关系。它不证明任何 residue 是实验表位、受体界面或中和位点。

### 10. 是否存在问题
UniProt 与冻结项目参考的差异已完整保留在映射中；项目分析继续使用 1–181 坐标。若 `pdb_structure_linkage_allowed=false`，所有依赖该 PDB 的结构联动必须立即停止。

### 11. 下一步
在 PDB/reference 门槛通过后，Agent B 继续尝试官方 standalone BepiPred-3.0 `vt_pred` 推理；结构 Agent 只可在该门槛通过时继续。

## Runtime provenance
- Python {platform.python_version()}
- Biopython {Bio.__version__}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report, encoding="utf-8")


def run_sequence_qc(
    reference_path: Path,
    pdb_path: Path,
    uniprot_path: Path,
    results_dir: Path,
    report_path: Path,
) -> dict:
    reference = read_fasta(reference_path)
    uniprot = read_fasta(uniprot_path)
    pdb = extract_chain(pdb_path, expected_chain="A")

    pdb_rows, pdb_summary = align_to_reference(reference.sequence, pdb.sequence, "pdb")
    for row in pdb_rows:
        subject_position = row["pdb_position"]
        if subject_position == "":
            row.update({"pdb_chain": "", "pdb_residue_id": "", "pdb_insertion_code": ""})
        else:
            index = int(subject_position) - 1
            row.update(
                {
                    "pdb_chain": pdb.chain_id,
                    "pdb_residue_id": pdb.pdb_residue_ids[index],
                    "pdb_insertion_code": pdb.pdb_insertion_codes[index],
                }
            )

    uniprot_rows, uniprot_summary = align_to_reference(
        reference.sequence, uniprot.sequence, "uniprot"
    )

    checks = {
        "reference_length_is_181": len(reference.sequence) == 181,
        "pdb_only_expected_chain": pdb.all_chain_ids == ["A"],
        "pdb_sequence_length_is_181": len(pdb.sequence) == 181,
        "pdb_matches_reference": pdb.sequence == reference.sequence,
        "residue_ids_map_to_project_positions": (
            pdb.project_positions == list(range(1, 182))
            and len(pdb.pdb_residue_ids) == 181
            and not any(pdb.pdb_insertion_codes)
        ),
        "no_missing_residues_or_backbone_atoms": (
            not pdb.missing_project_positions and not pdb.missing_backbone_atoms
        ),
    }
    qc = {
        "coordinate_system": "project_1_181",
        "reference_identifier": reference.identifier,
        "reference_length": len(reference.sequence),
        "pdb_chain_ids": pdb.all_chain_ids,
        "pdb_sequence_length": len(pdb.sequence),
        "pdb_residue_number_offsets": sorted(
            {
                residue_number - project_position
                for project_position, residue_number in zip(
                    pdb.project_positions, pdb.pdb_residue_numbers
                )
            }
        ),
        "uniprot_identifier": uniprot.identifier,
        "uniprot_length": len(uniprot.sequence),
        "checks": checks,
        "pdb_alignment": pdb_summary,
        "uniprot_alignment": uniprot_summary,
        "pdb_missing_project_positions": pdb.missing_project_positions,
        "pdb_missing_backbone_atoms": pdb.missing_backbone_atoms,
        "pdb_structure_linkage_allowed": all(checks.values()),
        "software": {"python": platform.python_version(), "biopython": Bio.__version__},
    }

    results_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(results_dir / "sequence_mapping.csv", pdb_rows)
    _write_fasta(
        results_dir / "pdb_sequence.fasta",
        "AF-Q925S4-F1-model_v6_chain_A project_mapping=1-181",
        pdb.sequence,
    )
    _write_csv(results_dir / "uniprot_mapping.csv", uniprot_rows)
    (results_dir / "sequence_qc.json").write_text(
        json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _write_report(report_path, qc, pdb_summary, uniprot_summary)
    return qc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--uniprot", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    qc = run_sequence_qc(
        args.reference, args.pdb, args.uniprot, args.results_dir, args.report
    )
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    return 0 if qc["pdb_structure_linkage_allowed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
