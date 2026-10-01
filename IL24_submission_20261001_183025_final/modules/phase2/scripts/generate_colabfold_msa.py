#!/usr/bin/env python3
"""Generate paired/unpaired ColabFold MSAs as Boltz CSV inputs.

This runs on the local workstation, where api.colabfold.com has a valid TLS
chain, and writes auditable CSV files that can be uploaded with the project.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "boltz" / "src"))

from boltz.data.const import max_msa_seqs, max_paired_seqs  # noqa: E402
from boltz.data.msa.mmseqs2 import run_mmseqs2  # noqa: E402


DATA = ROOT / "data"
OUT = DATA / "msa"
BOLTZ_CLI_MAX_MSA_SEQS = 8192


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name and line:
            records[name] += line.upper()
    return records


def combine_to_boltz_csv(paired_a3m: str, unpaired_a3m: str, output: Path) -> int:
    paired = paired_a3m.strip().splitlines()[1::2]
    paired = paired[:max_paired_seqs]
    keys = [idx for idx, seq in enumerate(paired) if seq != "-" * len(seq)]
    paired = [seq for seq in paired if seq != "-" * len(seq)]

    unpaired = unpaired_a3m.strip().splitlines()[1::2]
    unpaired = unpaired[: max_msa_seqs - len(paired)]
    if paired:
        unpaired = unpaired[1:]

    sequences = paired + unpaired
    keys.extend([-1] * len(unpaired))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["key", "sequence"])
        writer.writerows(zip(keys, sequences))
    return len(sequences)


def generate_target(name: str, antigen: str, vh: str, vl: str) -> list[dict[str, object]]:
    sequences = [antigen, vh, vl]
    work = OUT / name / "colabfold_raw"
    work.mkdir(parents=True, exist_ok=True)
    paired = run_mmseqs2(
        sequences,
        prefix=str(work / "paired"),
        use_env=True,
        use_pairing=True,
        pairing_strategy="greedy",
    )
    unpaired = run_mmseqs2(
        sequences,
        prefix=str(work / "unpaired"),
        use_env=True,
        use_pairing=False,
        pairing_strategy="greedy",
    )

    rows: list[dict[str, object]] = []
    for chain, paired_text, unpaired_text in zip("AHL", paired, unpaired):
        path = OUT / name / f"{chain}.csv"
        depth = combine_to_boltz_csv(paired_text, unpaired_text, path)
        rows.append(
            {
                "target": name,
                "chain": chain,
                "sequence_length": len(sequences["AHL".index(chain)]),
                "msa_depth": depth,
                "boltz_cli_max_msa_seqs": BOLTZ_CLI_MAX_MSA_SEQS,
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    return rows


def main() -> None:
    antibody = read_fasta(DATA / "VH_VL_clean.fasta")
    vh = next(seq for key, seq in antibody.items() if "VH" in key)
    vl = next(seq for key, seq in antibody.items() if "VL" in key)
    native = next(iter(read_fasta(DATA / "IL24_native_27_181.fasta").values()))
    immunogen = next(iter(read_fasta(DATA / "IL24_immunogen_27_160.fasta").values()))

    rows = generate_target("native", native, vh, vl)
    rows.extend(generate_target("immunogen", immunogen, vh, vl))
    manifest = OUT / "msa_manifest.json"
    manifest.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        print(
            f"{row['target']} chain {row['chain']}: depth={row['msa_depth']} "
            f"sha256={row['sha256']}"
        )
    print(f"Wrote manifest: {manifest}")


if __name__ == "__main__":
    main()
