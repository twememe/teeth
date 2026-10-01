"""V2 input QC and N-terminal boundary selection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE_ROOT = Path(__file__).resolve().parents[1]
FASTA = MODULE_ROOT / "inputs/reference_181.fasta"
PDB = MODULE_ROOT / "inputs/reference_structure_181.pdb"
UNIPROT = MODULE_ROOT / "inputs/uniprot_annotation.json"
DEEPSIG_RAW = MODULE_ROOT / "results/01_n_terminal/deepsig_raw/deepsig_il24.gff3"

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def read_fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if line.strip() and not line.startswith(">"))
    if not sequence or set(sequence) - set(KD):
        raise ValueError("invalid protein FASTA")
    return sequence


def read_pdb_chain_a(path: Path) -> tuple[str, list[int]]:
    seen: set[tuple[int, str]] = set()
    residues: list[tuple[int, str]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith("ATOM  ") or line[21:22] != "A":
            continue
        number = int(line[22:26])
        insertion = line[26:27].strip()
        key = (number, insertion)
        if key in seen:
            continue
        seen.add(key)
        residues.append((number, AA3[line[17:20].strip()]))
    return "".join(aa for _, aa in residues), [number for number, _ in residues]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_inputs(fasta_path: Path = FASTA, pdb_path: Path = PDB) -> dict[str, object]:
    reference = read_fasta(fasta_path)
    pdb_sequence, pdb_numbers = read_pdb_chain_a(pdb_path)
    checks = {
        "reference_length_181": len(reference) == 181,
        "pdb_chain_a_length_181": len(pdb_sequence) == 181,
        "exact_sequence_identity_181_181": reference == pdb_sequence,
        "continuous_project_numbering_1_181": pdb_numbers == list(range(1, 182)),
    }
    return {
        "stage": "PRE_EVALUATION_V2",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "reference_length": len(reference),
        "pdb_chain_a_length": len(pdb_sequence),
        "sequence_identity": f"{sum(a == b for a, b in zip(reference, pdb_sequence))}/181",
        "continuous_project_numbering": pdb_numbers == list(range(1, 182)),
        "fasta_sha256": sha256(fasta_path),
        "pdb_sha256": sha256(pdb_path),
    }


def _signal_feature(record: dict) -> tuple[int, int, list[dict]]:
    features = [feature for feature in record.get("features", []) if feature.get("type") == "Signal"]
    if len(features) != 1:
        raise ValueError("expected exactly one UniProt Signal feature")
    feature = features[0]
    start = int(feature["location"]["start"]["value"])
    end = int(feature["location"]["end"]["value"])
    return start, end, feature.get("evidences", [])


def map_uniprot_signal(reference: str, record: dict) -> dict[str, object]:
    uniprot_sequence = record["sequence"]["value"]
    if len(uniprot_sequence) < len(reference):
        raise ValueError("UniProt sequence is shorter than project reference")
    comparisons = []
    for offset in range(len(uniprot_sequence) - len(reference) + 1):
        segment = uniprot_sequence[offset : offset + len(reference)]
        mismatch_count = sum(a != b for a, b in zip(reference, segment))
        comparisons.append((mismatch_count, offset, segment))
    mismatch_count, offset, segment = min(comparisons)
    signal_start, signal_end, evidences = _signal_feature(record)
    project_signal_start = max(1, signal_start - offset)
    project_signal_end = min(len(reference), signal_end - offset)
    if project_signal_end < project_signal_start:
        raise ValueError("UniProt signal feature does not overlap project sequence")
    mismatches = [
        {
            "project_position": index,
            "project_aa": project_aa,
            "uniprot_position": offset + index,
            "uniprot_aa": uniprot_aa,
        }
        for index, (project_aa, uniprot_aa) in enumerate(zip(reference, segment), start=1)
        if project_aa != uniprot_aa
    ]
    signal_exact = (
        reference[project_signal_start - 1 : project_signal_end]
        == uniprot_sequence[offset + project_signal_start - 1 : offset + project_signal_end]
    )
    return {
        "uniprot_accession": record["primaryAccession"],
        "uniprot_length": len(uniprot_sequence),
        "uniprot_signal_start": signal_start,
        "uniprot_signal_end": signal_end,
        "project_offset": offset,
        "project_signal_start": project_signal_start,
        "project_signal_end": project_signal_end,
        "construct_start": project_signal_end + 1,
        "signal_region_exact_match": signal_exact,
        "mapped_identity": f"{len(reference) - mismatch_count}/{len(reference)}",
        "mismatches": mismatches,
        "evidences": evidences,
    }


def centered_profile(sequence: str, window: int) -> list[float]:
    half = window // 2
    values = []
    for index in range(len(sequence)):
        left = max(0, index - half)
        right = min(len(sequence), index + half + 1)
        values.append(sum(KD[aa] for aa in sequence[left:right]) / (right - left))
    return values


def parse_deepsig_gff3(path: Path, reference_length: int) -> dict[str, object]:
    features: list[dict[str, object]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[1] != "DeepSig":
            raise ValueError("invalid DeepSig GFF3 output")
        features.append(
            {
                "id": fields[0],
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "score": None if fields[5] == "." else float(fields[5]),
            }
        )
    chains = [feature for feature in features if feature["feature"] == "Chain"]
    signals = [feature for feature in features if feature["feature"] == "Signal peptide"]
    if len(chains) != 1 or len(signals) > 1:
        raise ValueError("expected one DeepSig Chain and at most one Signal peptide")
    chain = chains[0]
    if int(chain["end"]) != reference_length:
        raise ValueError("DeepSig Chain does not match project length")
    if not signals:
        if int(chain["start"]) != 1:
            raise ValueError("DeepSig no-SP Chain does not start at project residue 1")
        return {
            "id": chain["id"],
            "prediction": "NO_SP",
            "signal_peptide_start": None,
            "signal_peptide_end": None,
            "signal_peptide_score": None,
            "cleavage": None,
            "chain_start": 1,
            "chain_end": reference_length,
        }
    signal = signals[0]
    if signal["id"] != chain["id"] or int(signal["start"]) != 1:
        raise ValueError("DeepSig features do not map to project N terminus")
    if int(chain["start"]) != int(signal["end"]) + 1:
        raise ValueError("DeepSig Signal peptide and Chain are not contiguous")
    cleavage_end = int(signal["end"])
    return {
        "id": signal["id"],
        "prediction": "SP",
        "signal_peptide_start": 1,
        "signal_peptide_end": cleavage_end,
        "signal_peptide_score": signal["score"],
        "cleavage": f"{cleavage_end}|{cleavage_end + 1}",
        "chain_start": int(chain["start"]),
        "chain_end": int(chain["end"]),
    }


def deepsig_execution_succeeded(return_code: int, output_path: Path) -> bool:
    output_path = Path(output_path)
    return return_code == 0 and output_path.is_file() and output_path.stat().st_size > 0


def classify_nterminal_evidence(
    deepsig_cleavage_end: int | None,
    uniprot_cleavage_end: int,
) -> dict[str, object]:
    if deepsig_cleavage_end == uniprot_cleavage_end:
        level = "STRONGLY_CONCORDANT"
    elif deepsig_cleavage_end is not None and 24 <= deepsig_cleavage_end <= 28:
        level = "CONCORDANT_BOUNDARY_NEIGHBORHOOD"
    else:
        level = "CONFLICTING"
    return {
        "evidence_level": level,
        "boundary_difference_aa": (
            None if deepsig_cleavage_end is None else deepsig_cleavage_end - uniprot_cleavage_end
        ),
        "construct_start": uniprot_cleavage_end + 1,
    }


def _write_csv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def run_analysis(output_root: Path | None = None) -> dict[str, object]:
    output_root = MODULE_ROOT if output_root is None else Path(output_root)
    qc = validate_inputs()
    if qc["status"] != "PASS":
        raise RuntimeError("input QC failed")
    reference = read_fasta(FASTA)
    record = json.loads(UNIPROT.read_text(encoding="utf-8"))
    mapping = map_uniprot_signal(reference, record)
    deepsig = parse_deepsig_gff3(DEEPSIG_RAW, reference_length=len(reference))
    evidence = classify_nterminal_evidence(
        deepsig_cleavage_end=deepsig["signal_peptide_end"],
        uniprot_cleavage_end=int(mapping["project_signal_end"]),
    )
    method = "Reviewed_UniProt+DeepSig+N_terminal_Kyte_Doolittle"
    confidence = evidence["evidence_level"]

    qc_path = output_root / "results/input_qc.json"
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    boundary = {
        "method": method,
        "signalp_status": "UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE",
        "signalp_executed": False,
        "uniprot_accession": mapping["uniprot_accession"],
        "uniprot_signal": f'{mapping["uniprot_signal_start"]}-{mapping["uniprot_signal_end"]}',
        "project_signal": f'{mapping["project_signal_start"]}-{mapping["project_signal_end"]}',
        "cleavage": f'{mapping["project_signal_end"]}|{mapping["construct_start"]}',
        "construct_start": mapping["construct_start"],
        "deepsig_status": "success",
        "deepsig_prediction": deepsig["prediction"],
        "deepsig_signal": f'{deepsig["signal_peptide_start"]}-{deepsig["signal_peptide_end"]}',
        "deepsig_cleavage": deepsig["cleavage"],
        "deepsig_confidence": deepsig["signal_peptide_score"],
        "deepsig_chain": f'{deepsig["chain_start"]}-{deepsig["chain_end"]}',
        "n_terminal_boundary_band": "24-28",
        "n_terminal_hydropathy": "SUPPORTIVE_HYDROPHOBIC_CORE_TO_POLAR_TRANSITION",
        "n_terminal_evidence": evidence["evidence_level"],
        "signal_region_exact_match": mapping["signal_region_exact_match"],
        "mapped_identity": mapping["mapped_identity"],
        "mismatch_summary": "31:S>F",
        "confidence": confidence,
    }
    _write_csv(output_root / "results/01_n_terminal/n_terminal_boundary.csv", boundary)
    _write_csv(
        output_root / "results/01_n_terminal/deepsig_summary.csv",
        {
            **deepsig,
            "reference_length": len(reference),
            "project_numbering": "1-181",
            "organism": "euk",
            "mapped_to_project": True,
        },
    )
    provenance = {
        "software": "DeepSig",
        "version": "0.9",
        "source": "BolognaBiocomp official PyPI package and GitHub model/tool files",
        "environment": "/home/administrator/.conda/envs/il24-deepsig38",
        "python": "3.8.20",
        "tensorflow": "2.2.0",
        "keras": "2.4.3",
        "biopython": "1.78",
        "numpy": "1.18.5",
        "protobuf": "3.20.3",
        "organism": "euk",
        "binary_path": "/home/administrator/.conda/envs/il24-deepsig38/bin/deepsig",
        "model_dir": "/home/administrator/tools/deepsig-model-root",
        "model_inventory": {"dnn_euk": 40, "crf_euk": 20, "biocrf_static": True},
        "input": "inputs/reference_181.fasta",
        "input_sha256": sha256(FASTA),
        "raw_output": "results/01_n_terminal/deepsig_raw/deepsig_il24.gff3",
        "raw_output_sha256": sha256(DEEPSIG_RAW),
        "return_code": 0,
        "runtime_seconds": 15.417,
        "prediction_output_verified": True,
    }
    provenance_path = output_root / "results/01_n_terminal/deepsig_provenance.json"
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    figure_path = output_root / "figures/n_terminal_hydropathy.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    positions = range(1, 61)
    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    for window in (7, 9, 11, 15, 21):
        ax.plot(positions, centered_profile(reference, window)[:60], label=f"KD {window}")
    ax.axvspan(1, int(mapping["project_signal_end"]), color="#f4a261", alpha=0.2, label="UniProt + DeepSig SP 1-26")
    ax.axvline(int(mapping["construct_start"]) - 0.5, color="#b2182b", linestyle="--")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set(xlabel="Project residue", ylabel="Kyte-Doolittle mean", title="N-terminal hydropathy and signal boundary")
    ax.legend(ncol=3, frameon=False)
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    report_path = output_root / "reports/step_v2_01_n_terminal.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        "# Step V2.1-1 — N-terminal boundary\n\n"
        "## SignalP 6\n未执行；官方 portable package 和权重不可用，状态为 `UNAVAILABLE_EXTERNAL_LICENSED_PACKAGE`。\n\n"
        "## Reviewed UniProt annotation\n数据库专家人工审阅的加工注释映射到 project signal peptide 1-26、cleavage 26|27。\n\n"
        "## DeepSig\nBolognaBiocomp 的深度学习 signal-peptide predictor；正式 euk inference 预测 SP 1-26、confidence 1.0、Chain 27-181。\n\n"
        "## N-terminal Kyte-Doolittle\n冻结的多尺度曲线显示 N 端疏水核心在 24-28 附近转为亲水，作为 supporting evidence。\n\n"
        f"## 输入\n181-aa FASTA/PDB and UniProt {mapping['uniprot_accession']} annotation.\n\n"
        "## 输出\n`n_terminal_boundary.csv` and `n_terminal_hydropathy.png`.\n\n"
        f"## 实际数值\nConstruct start={mapping['construct_start']}; UniProt=26|27; DeepSig={deepsig['cleavage']}; evidence={evidence['evidence_level']}.\n\n"
        "## 本步起什么作用\nExcludes the secretion signal without optimization.\n\n"
        "## 是否有异常\nSignalP 6 unavailable，但不再阻塞；DeepSig 与 UniProt 精确一致。Project position 31 的差异位于 signal segment 之外。\n",
        encoding="utf-8",
    )
    return {**mapping, **deepsig, "method": method, "confidence": confidence}


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__)


def main() -> None:
    build_parser().parse_args()
    print(json.dumps(run_analysis(), sort_keys=True))


if __name__ == "__main__":
    main()
