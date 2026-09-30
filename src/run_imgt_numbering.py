#!/usr/bin/env python3
"""Run ANARCI IMGT numbering and write auditable residue/CDR tables."""
from __future__ import annotations

import csv
import os
from pathlib import Path

from anarci import anarci

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data" / "VH_VL_clean.fasta"
OUTPUT_ROOT = Path(os.environ.get("PHASE2_OUTPUT_ROOT", ROOT / "results")).resolve()
HMMER = Path(os.environ.get("PHASE2_HMMER_BIN", ""))
OUT = OUTPUT_ROOT / "10_numbering"
RESIDUE_OUT = OUT / "antibody_imgt_numbering.csv"
CDR_OUT = OUT / "cdr_annotation.csv"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = None
    parts: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line.upper())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def region_for_imgt(position: int) -> str:
    if position <= 26:
        return "FR1"
    if position <= 38:
        return "CDR1"
    if position <= 55:
        return "FR2"
    if position <= 65:
        return "CDR2"
    if position <= 104:
        return "FR3"
    if position <= 117:
        return "CDR3"
    return "FR4"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = read_fasta(FASTA)
    numbered, details, _ = anarci(
        records,
        scheme="imgt",
        hmmerpath=str(HMMER),
        ncpu=2,
        assign_germline=False,
    )

    residue_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for (name, sequence), domains, domain_details in zip(records, numbered, details):
        if not domains:
            raise RuntimeError(f"ANARCI found no antibody domain for {name}")
        if len(domains) != 1:
            raise RuntimeError(f"Expected one domain for {name}, found {len(domains)}")

        numbering, domain_start, domain_end = domains[0]
        chain_type = domain_details[0]["chain_type"]
        chain = "VH" if chain_type == "H" else "VL"
        sequence_index = domain_start + 1
        by_region: dict[str, list[tuple[int, str, str, int]]] = {}

        for (position, insertion), aa in numbering:
            if aa == "-":
                continue
            region = region_for_imgt(position)
            by_region.setdefault(region, []).append((position, insertion.strip(), aa, sequence_index))
            residue_rows.append({
                "antibody_id": "IA6-13-8",
                "chain": chain,
                "sequence_index": sequence_index,
                "imgt_position": position,
                "imgt_insertion": insertion.strip(),
                "amino_acid": aa,
                "region": region,
            })
            if sequence[sequence_index - 1] != aa:
                raise RuntimeError(
                    f"Sequence mismatch for {name} at index {sequence_index}: "
                    f"input={sequence[sequence_index - 1]} ANARCI={aa}"
                )
            sequence_index += 1

        if sequence_index - 1 != domain_end + 1:
            raise RuntimeError(f"ANARCI domain end mismatch for {name}")

        for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
            residues = by_region.get(region, [])
            if not residues:
                continue
            first = residues[0]
            last = residues[-1]
            summary_rows.append({
                "antibody_id": "IA6-13-8",
                "chain": chain,
                "region": region,
                "sequence_start": first[3],
                "sequence_end": last[3],
                "imgt_start": f"{first[0]}{first[1]}",
                "imgt_end": f"{last[0]}{last[1]}",
                "length": len(residues),
                "sequence": "".join(item[2] for item in residues),
            })

    with RESIDUE_OUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(residue_rows[0]))
        writer.writeheader()
        writer.writerows(residue_rows)

    with CDR_OUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {RESIDUE_OUT} ({len(residue_rows)} residues)")
    print(f"Wrote {CDR_OUT} ({len(summary_rows)} regions)")


if __name__ == "__main__":
    main()
