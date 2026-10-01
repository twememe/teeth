#!/usr/bin/env python
"""Mammalian IL-24 ortholog selection and residue conservation analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import subprocess
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MAMMALIAN_ORDERS = (
    "Afrosoricida",
    "Artiodactyla",
    "Carnivora",
    "Chiroptera",
    "Cingulata",
    "Dasyuromorphia",
    "Didelphimorphia",
    "Diprotodontia",
    "Eulipotyphla",
    "Hyracoidea",
    "Lagomorpha",
    "Macroscelidea",
    "Monotremata",
    "Perissodactyla",
    "Pholidota",
    "Pilosa",
    "Primates",
    "Proboscidea",
    "Rodentia",
    "Scandentia",
    "Sirenia",
    "Tubulidentata",
)


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = None
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records[identifier] = "".join(sequence).upper()
            identifier = line[1:].split()[0]
            if identifier in records:
                raise ValueError(f"Duplicate FASTA identifier: {identifier}")
            sequence = []
        elif identifier is None:
            raise ValueError(f"Sequence precedes FASTA header in {path}")
        else:
            sequence.append(line)
    if identifier is not None:
        records[identifier] = "".join(sequence).upper()
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records


def load_reference(path: Path) -> tuple[str, str]:
    records = parse_fasta(path)
    if len(records) != 1:
        raise ValueError(f"Expected one frozen reference record in {path}")
    identifier, sequence = next(iter(records.items()))
    if set(sequence) - AMINO_ACIDS:
        raise ValueError("Frozen reference contains non-canonical amino acids")
    return identifier, sequence


def load_uniprot_json(path: Path) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload.get("results", [])
    if not records:
        raise ValueError(f"No UniProt results in {path}")
    return records


def _protein_name(description: dict) -> str:
    recommended = description.get("recommendedName")
    if recommended:
        return recommended["fullName"]["value"]
    submitted = description.get("submissionNames") or []
    if submitted:
        return submitted[0]["fullName"]["value"]
    return ""


def _mammalian_order(lineage: list[str]) -> str:
    for order in MAMMALIAN_ORDERS:
        if order in lineage:
            return order
    return "unclassified"


def _compact_record(raw: dict) -> dict:
    description = raw.get("proteinDescription") or {}
    genes = raw.get("genes") or []
    gene_name = ""
    if genes and genes[0].get("geneName"):
        gene_name = genes[0]["geneName"]["value"]
    organism = raw["organism"]
    sequence = raw["sequence"]["value"].upper()
    entry_type = raw.get("entryType", "")
    accession = raw["primaryAccession"]
    return {
        "accession": accession,
        "uniprot_id": raw.get("uniProtkbId", ""),
        "species": organism["scientificName"],
        "taxid": int(organism["taxonId"]),
        "mammalian_order": _mammalian_order(organism.get("lineage", [])),
        "reviewed": entry_type.startswith("UniProtKB reviewed"),
        "annotation_score": float(raw.get("annotationScore") or 0.0),
        "protein_name": _protein_name(description),
        "protein_flag": str(description.get("flag") or ""),
        "gene_name": gene_name,
        "length": len(sequence),
        "sequence": sequence,
        "source_url": f"https://rest.uniprot.org/uniprotkb/{accession}",
    }


def _eligibility_reasons(record: dict) -> list[str]:
    reasons: list[str] = []
    name = record["protein_name"].lower()
    flag = record["protein_flag"].lower()
    if "low quality" in name:
        reasons.append("explicit_low_quality_name")
    if "isoform" in name:
        reasons.append("isoform_name")
    if "fragment" in name or flag == "fragment":
        reasons.append("explicit_fragment")
    if record["length"] < 140:
        reasons.append("length_below_140")
    if record["length"] > 280:
        reasons.append("length_above_280")
    if set(record["sequence"]) - AMINO_ACIDS:
        reasons.append("noncanonical_amino_acid")
    return reasons


def _quality_key(record: dict) -> tuple:
    normalized_name = record["protein_name"].lower().replace("-", " ")
    canonical_name = normalized_name == "interleukin 24"
    exact_gene = record["gene_name"].lower() == "il24"
    return (
        -int(record["reviewed"]),
        -int(exact_gene),
        -int(canonical_name),
        -record["annotation_score"],
        abs(record["length"] - 200),
        record["accession"],
    )


def select_representatives(
    raw_records: list[dict], target_count: int = 25
) -> tuple[list[dict], list[dict]]:
    if not 15 <= target_count <= 30:
        raise ValueError("target_count must be between 15 and 30")
    compact = [_compact_record(raw) for raw in raw_records]
    audit_by_accession: dict[str, dict] = {}
    eligible_by_taxid: dict[int, list[dict]] = defaultdict(list)

    for record in compact:
        reasons = _eligibility_reasons(record)
        if reasons:
            audit_by_accession[record["accession"]] = {
                **record,
                "selection_status": "excluded",
                "selection_reason": ";".join(reasons),
            }
        else:
            eligible_by_taxid[record["taxid"]].append(record)

    species_representatives: list[dict] = []
    for taxid in sorted(eligible_by_taxid):
        ranked = sorted(eligible_by_taxid[taxid], key=_quality_key)
        chosen = ranked[0]
        species_representatives.append(chosen)
        for duplicate in ranked[1:]:
            audit_by_accession[duplicate["accession"]] = {
                **duplicate,
                "selection_status": "excluded",
                "selection_reason": f"same_taxid_lower_priority_than_{chosen['accession']}",
            }

    buckets: dict[str, list[dict]] = defaultdict(list)
    for record in species_representatives:
        buckets[record["mammalian_order"]].append(record)
    for order in buckets:
        buckets[order].sort(key=_quality_key)

    selected: list[dict] = []
    depth = 0
    orders = sorted(buckets)
    while len(selected) < target_count:
        added = False
        for order in orders:
            if depth < len(buckets[order]):
                selected.append(buckets[order][depth])
                added = True
                if len(selected) == target_count:
                    break
        if not added:
            break
        depth += 1
    if len(selected) < target_count:
        raise ValueError(
            f"Only {len(selected)} eligible species representatives; requested {target_count}"
        )

    selected_accessions = {record["accession"] for record in selected}
    for record in species_representatives:
        included = record["accession"] in selected_accessions
        audit_by_accession[record["accession"]] = {
            **record,
            "selection_status": "included" if included else "excluded",
            "selection_reason": (
                "order_round_robin_quality_representative"
                if included
                else f"eligible_not_selected_target_{target_count}"
            ),
        }
    selected = [
        {
            **record,
            "selection_reason": "order_round_robin_quality_representative",
        }
        for record in selected
    ]
    audit = [audit_by_accession[record["accession"]] for record in compact]
    return selected, audit


def _write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, description, sequence in records:
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def write_ortholog_fasta(path: Path, selected: list[dict]) -> None:
    records = [
        (
            record["accession"],
            f"taxid={record['taxid']} species={record['species'].replace(' ', '_')}",
            record["sequence"],
        )
        for record in selected
    ]
    _write_fasta(path, records)


def build_alignment_input(
    path: Path,
    reference_id: str,
    reference_sequence: str,
    selected: list[dict],
) -> None:
    records = [(reference_id, "coordinate_system=project_1_181", reference_sequence)]
    records.extend(
        (
            record["accession"],
            f"taxid={record['taxid']} species={record['species'].replace(' ', '_')}",
            record["sequence"],
        )
        for record in selected
    )
    _write_fasta(path, records)


def run_mafft(
    input_path: Path,
    output_path: Path,
    mafft_executable: str = "mafft",
    log_path: Path | None = None,
) -> None:
    command = [mafft_executable, "--auto", "--thread", "1", str(input_path)]
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
    )
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            f"command={' '.join(command)}\nreturncode={completed.returncode}\n\n{completed.stderr}",
            encoding="utf-8",
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"MAFFT exited {completed.returncode}: {completed.stderr.strip()}"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(completed.stdout, encoding="utf-8")


def _reference_columns(alignment: dict[str, str], reference_id: str) -> list[int]:
    if reference_id not in alignment:
        raise ValueError(f"Reference {reference_id!r} absent from alignment")
    lengths = {len(sequence) for sequence in alignment.values()}
    if len(lengths) != 1:
        raise ValueError("Aligned FASTA records do not share one length")
    return [index for index, residue in enumerate(alignment[reference_id]) if residue != "-"]


def _column_conservation(residues: list[str]) -> tuple[float, float, float, int]:
    non_gaps = [residue for residue in residues if residue != "-"]
    gap_fraction = 1.0 - len(non_gaps) / len(residues)
    if not non_gaps:
        return 0.0, 1.0, gap_fraction, 0
    counts = Counter(non_gaps)
    entropy = -sum(
        (count / len(non_gaps)) * math.log(count / len(non_gaps))
        for count in counts.values()
    )
    normalized_entropy = entropy / math.log(20)
    conservation = max(0.0, min(1.0, 1.0 - normalized_entropy))
    return conservation, entropy, gap_fraction, len(non_gaps)


def compute_conservation(
    alignment: dict[str, str], reference_id: str
) -> list[dict]:
    reference = alignment[reference_id]
    orthologs = [sequence for key, sequence in alignment.items() if key != reference_id]
    if not orthologs:
        raise ValueError("Alignment contains no ortholog sequences")
    columns = _reference_columns(alignment, reference_id)
    rows: list[dict] = []
    for position, column in enumerate(columns, start=1):
        conservation, entropy, gap_fraction, non_gap_count = _column_conservation(
            [sequence[column] for sequence in orthologs]
        )
        rows.append(
            {
                "position": position,
                "AA": reference[column],
                "alignment_column": column + 1,
                "conservation_score": conservation,
                "shannon_entropy": entropy,
                "normalized_entropy": 1.0 - conservation,
                "gap_fraction": gap_fraction,
                "non_gap_count": non_gap_count,
                "ortholog_count": len(orthologs),
            }
        )
    return rows


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    index = (len(ordered) - 1) * fraction
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    weight = index - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def bootstrap_conservation(
    alignment: dict[str, str],
    reference_id: str,
    replicates: int = 100,
    seed: int = 2406,
) -> list[dict]:
    if replicates != 100:
        raise ValueError("This protocol fixes conservation bootstrap replicates at 100")
    orthologs = [sequence for key, sequence in alignment.items() if key != reference_id]
    columns = _reference_columns(alignment, reference_id)
    rng = random.Random(seed)
    replicate_scores = [[] for _ in columns]
    for _ in range(replicates):
        sampled = [orthologs[rng.randrange(len(orthologs))] for _ in orthologs]
        for row_index, column in enumerate(columns):
            score, _, _, _ = _column_conservation(
                [sequence[column] for sequence in sampled]
            )
            replicate_scores[row_index].append(score)
    return [
        {
            "position": position,
            "bootstrap_median": statistics.median(scores),
            "bootstrap_ci_lower": _percentile(scores, 0.025),
            "bootstrap_ci_upper": _percentile(scores, 0.975),
            "bootstrap_replicates": replicates,
            "bootstrap_seed": seed,
        }
        for position, scores in enumerate(replicate_scores, start=1)
    ]


def _write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_profile(path: Path, rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = [row["position"] for row in rows]
    score = [row["conservation_score"] for row in rows]
    lower = [row["bootstrap_ci_lower"] for row in rows]
    upper = [row["bootstrap_ci_upper"] for row in rows]
    gaps = [row["gap_fraction"] for row in rows]
    figure, axis = plt.subplots(figsize=(11, 4.8))
    axis.plot(x, score, color="#2166ac", linewidth=1.4, label="Normalized conservation")
    axis.fill_between(x, lower, upper, color="#67a9cf", alpha=0.25, label="95% bootstrap interval")
    axis.plot(x, gaps, color="#b2182b", linewidth=1.0, alpha=0.8, label="Gap fraction")
    axis.set(xlabel="Project position (1–181)", ylabel="Score / fraction", xlim=(1, 181), ylim=(0, 1.02))
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, ncol=3, fontsize=8)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200)
    plt.close(figure)


def write_report(
    path: Path,
    selected: list[dict],
    audit: list[dict],
    rows: list[dict],
    mafft_version: str,
) -> None:
    scores = [row["conservation_score"] for row in rows]
    gaps = [row["gap_fraction"] for row in rows]
    excluded = Counter(
        reason
        for row in audit
        if row["selection_status"] == "excluded"
        for reason in row["selection_reason"].split(";")
    )
    orders = sorted({record["mammalian_order"] for record in selected})
    report = f"""# Step 6 — Evolutionary Conservation

## Step 6 完成

### 1. 本步做了什么
从 UniProt 官方 REST API 的 mammalian IL-24 查询结果中，按预注册质量规则选择 one species/one representative 的 25 条序列，与冻结的 181-aa project reference 一起用 MAFFT 比对，并逐位计算保守性、gap fraction 与 100 次固定种子 bootstrap 区间。

### 2. 使用的算法
UniProt REST API；MAFFT `{mafft_version}`；normalized Shannon entropy；100-replicate sequence bootstrap（seed=2406）。

### 3. 算法属于什么
UniProt 检索和 MAFFT 属于传统生物信息学；Shannon entropy 和 bootstrap 属于统计方法。没有训练或微调 AI 模型。

### 4. 算法原理
通俗地说，MAFFT 把不同物种的 IL-24 像多篇文章逐字排齐；某列越多物种使用同一种氨基酸，该位置越保守。技术上，去除 gap 后计算 `H=-sum(p*ln p)`，再转换为 `1-H/ln(20)`；gap 单独报告，不伪装成氨基酸。每次 bootstrap 对物种序列有放回抽样，并重新计算 181 个位置。

### 5. 为什么本项目需要它
它提供独立于序列表位、三维表位和表面可接近性的 evolutionary evidence，用于判断哪些 project residues 在哺乳动物中长期稳定。

### 6. 输入
- `inputs/reference_181.fasta`
- `results/06_conservation/raw/uniprot_il24_mammalia.json`

### 7. 输出
- `results/06_conservation/orthologs.fasta`
- `results/06_conservation/ortholog_metadata.csv`
- `results/06_conservation/ortholog_selection_audit.csv`
- `results/06_conservation/mafft_alignment.fasta`
- `results/06_conservation/residue_conservation.csv`
- `results/06_conservation/conservation_profile.png`
- `reports/step06_conservation.md`

### 8. 关键数值结果
- UniProt query records：{len(audit)}；selected species：{len(selected)}；reviewed selected：{sum(record['reviewed'] for record in selected)}。
- Mammalian orders ({len(orders)})：{', '.join(orders)}。
- Conservation min/mean/median/max：{min(scores):.6f} / {statistics.fmean(scores):.6f} / {statistics.median(scores):.6f} / {max(scores):.6f}。
- Gap fraction min/mean/max：{min(gaps):.6f} / {statistics.fmean(gaps):.6f} / {max(gaps):.6f}。
- Exclusion/audit reasons：{dict(excluded)}。

### 9. 当前结果怎样解释
高 conservation 表示该对齐位置在所选哺乳动物 IL-24 中氨基酸组成较一致；它不证明实验表位、受体界面或中和功能。高 gap fraction 位置需要降低解释置信度。

### 10. 是否存在问题
只有三条查询结果是 reviewed Swiss-Prot；其余代表来自透明筛选的 unreviewed entries。Bootstrap 反映当前 ortholog set 的抽样不确定性，不覆盖数据库注释错误或系统发育非独立性。

### 11. 下一步
Integration Agent 只读取冻结的 181-row conservation table，并在不访问历史区间的前提下与其他独立证据整合。
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--uniprot-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--target-count", type=int, default=25)
    parser.add_argument("--bootstrap-seed", type=int, default=2406)
    parser.add_argument("--retrieval-date", default=date.today().isoformat())
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reference_id, reference = load_reference(args.reference)
    if len(reference) != 181:
        raise ValueError(f"Frozen project reference length is {len(reference)}, expected 181")
    raw_records = load_uniprot_json(args.uniprot_json)
    selected, audit = select_representatives(raw_records, args.target_count)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    write_ortholog_fasta(output / "orthologs.fasta", selected)
    metadata_fields = [
        "accession", "species", "taxid", "reviewed", "length", "mammalian_order",
        "gene_name", "protein_name", "annotation_score", "source_url", "retrieval_date",
        "selection_reason",
    ]
    metadata = [{**record, "retrieval_date": args.retrieval_date} for record in selected]
    _write_csv(output / "ortholog_metadata.csv", metadata, metadata_fields)
    audit_fields = metadata_fields[:-2] + ["selection_status", "selection_reason"]
    audit_rows = [{**record, "retrieval_date": args.retrieval_date} for record in audit]
    _write_csv(output / "ortholog_selection_audit.csv", audit_rows, audit_fields)

    alignment_input = output / "raw" / "mafft_input_with_project_reference.fasta"
    mafft_log = output / "raw" / "mafft.log"
    build_alignment_input(alignment_input, reference_id, reference, selected)
    run_mafft(
        alignment_input,
        output / "mafft_alignment.fasta",
        mafft_executable=args.mafft,
        log_path=mafft_log,
    )
    alignment = parse_fasta(output / "mafft_alignment.fasta")
    table = compute_conservation(alignment, reference_id)
    if len(table) != 181 or "".join(row["AA"] for row in table) != reference:
        raise ValueError("MAFFT reference mapping did not produce project positions 1–181")
    bootstrap = bootstrap_conservation(
        alignment, reference_id, replicates=100, seed=args.bootstrap_seed
    )
    merged = [{**row, **bootstrap_row} for row, bootstrap_row in zip(table, bootstrap)]
    _write_csv(output / "residue_conservation.csv", merged)
    write_profile(output / "conservation_profile.png", merged)

    version = subprocess.run(
        [args.mafft, "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False
    ).stdout.strip()
    write_report(args.report, selected, audit, merged, version)
    qc = {
        "coordinate_system": "project_1_181",
        "query_record_count": len(raw_records),
        "selected_species_count": len(selected),
        "selected_taxids_unique": len({record["taxid"] for record in selected}) == len(selected),
        "selected_order_count": len({record["mammalian_order"] for record in selected}),
        "reference_length": len(reference),
        "alignment_sequence_count": len(alignment),
        "residue_rows": len(merged),
        "bootstrap_replicates": 100,
        "bootstrap_seed": args.bootstrap_seed,
        "mafft_version": version,
    }
    (output / "conservation_qc.json").write_text(
        json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
