#!/usr/bin/env python3
"""Freeze Phase 2 sequences and residue-number mappings from local inputs."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PDB = DATA / "AF-Q925S4-F1-model_v6.pdb"
ANTIBODY_FASTA = DATA / "VH_VL.fasta"

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_chain_sequence(path: Path, chain: str = "A") -> tuple[str, list[int]]:
    residues: list[tuple[int, str]] = []
    seen: set[tuple[int, str]] = set()
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            if line[21].strip() != chain:
                continue
            number = int(line[22:26])
            insertion = line[26].strip()
            key = (number, insertion)
            if key in seen:
                continue
            seen.add(key)
            residue = line[17:20].strip()
            if residue not in AA3_TO_1:
                raise ValueError(f"Unsupported residue {residue} at {chain}:{number}{insertion}")
            residues.append((number, AA3_TO_1[residue]))

    if not residues:
        raise ValueError(f"No CA residues found for chain {chain} in {path}")
    numbers = [number for number, _ in residues]
    return "".join(aa for _, aa in residues), numbers


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts).upper()
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    return records


def validate_protein(name: str, sequence: str) -> None:
    invalid = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid:
        raise ValueError(f"{name} contains invalid residues: {invalid}")


def write_fasta(path: Path, header: str, sequence: str) -> None:
    lines = [f">{header}"]
    lines.extend(sequence[i:i + 80] for i in range(0, len(sequence), 80))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    full, pdb_numbers = extract_chain_sequence(PDB, "A")
    if len(full) != 181:
        raise ValueError(f"Expected 181 IL-24 residues, found {len(full)}")
    if pdb_numbers != list(range(1, 182)):
        raise ValueError("PDB chain A numbering is not exactly 1..181")
    validate_protein("IL-24", full)

    antibodies = read_fasta(ANTIBODY_FASTA)
    if len(antibodies) != 2:
        raise ValueError(f"Expected two antibody FASTA records, found {len(antibodies)}")
    for name, sequence in antibodies.items():
        validate_protein(name, sequence)

    native = full[26:181]
    immunogen = full[26:160]
    if len(native) != 155 or len(immunogen) != 134:
        raise AssertionError("Unexpected mature/immunogen sequence lengths")

    full_path = DATA / "IL24_full_1_181.fasta"
    native_path = DATA / "IL24_native_27_181.fasta"
    immunogen_path = DATA / "IL24_immunogen_27_160.fasta"
    clean_ab_path = DATA / "VH_VL_clean.fasta"
    map_path = DATA / "il24_numbering_map.csv"
    manifest_path = DATA / "phase2_input_manifest_sha256.csv"

    write_fasta(full_path, "IL24_full_project_numbering_1-181", full)
    write_fasta(native_path, "IL24_native_project_27-181_model_1-155", native)
    write_fasta(immunogen_path, "IL24_immunogen_project_27-160_model_1-134", immunogen)

    clean_lines: list[str] = []
    for name, sequence in antibodies.items():
        clean_lines.append(f">{name}")
        clean_lines.extend(sequence[i:i + 80] for i in range(0, len(sequence), 80))
    clean_ab_path.write_text("\n".join(clean_lines) + "\n", encoding="ascii")

    with map_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "project_residue_number", "amino_acid", "native_27_181_index",
            "immunogen_27_160_index", "region"
        ])
        writer.writeheader()
        for project_number, aa in enumerate(full, 1):
            writer.writerow({
                "project_residue_number": project_number,
                "amino_acid": aa,
                "native_27_181_index": project_number - 26 if 27 <= project_number <= 181 else "",
                "immunogen_27_160_index": project_number - 26 if 27 <= project_number <= 160 else "",
                "region": (
                    "signal_peptide" if project_number <= 26 else
                    "immunogen" if project_number <= 160 else
                    "native_C_terminal_extension"
                ),
            })

    frozen = [PDB, ANTIBODY_FASTA, full_path, native_path, immunogen_path, clean_ab_path, map_path]
    with manifest_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
        writer.writeheader()
        for path in frozen:
            writer.writerow({
                "relative_path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })

    print(f"IL24 full: {len(full)} aa")
    print(f"IL24 native: {len(native)} aa")
    print(f"IL24 immunogen: {len(immunogen)} aa")
    for name, sequence in antibodies.items():
        print(f"{name}: {len(sequence)} aa")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
