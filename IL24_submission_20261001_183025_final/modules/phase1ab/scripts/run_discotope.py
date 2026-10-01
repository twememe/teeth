#!/usr/bin/env python3
"""Validated wrapper utilities for official DiscoTope-3.0 inference."""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


EXPECTED_LENGTH = 181
REQUIRED_OFFICIAL_COLUMNS = (
    "pdb",
    "res_id",
    "residue",
    "DiscoTope-3.0_score",
    "rsa",
    "pLDDTs",
    "length",
    "alphafold_struc_flag",
)

PHASE1B_ADDITIONAL_OFFICIAL_COLUMNS = (
    "chain",
    "calibrated_score",
    "epitope",
)


class MappingError(ValueError):
    """Raised when official output cannot map exactly to project coordinates."""


def load_reference(path: str | Path) -> str:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    sequence = sequence.upper()
    if len(sequence) != EXPECTED_LENGTH:
        raise MappingError(
            f"Frozen project reference must be {EXPECTED_LENGTH} residues; found {len(sequence)}"
        )
    return sequence


def official_repo_commit(repo_dir: str | Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def validate_project_input(
    pdb_path: str | Path, reference: str
) -> tuple[list[dict], dict]:
    structure = PDBParser(QUIET=True).get_structure("IL24_DT3", str(pdb_path))
    models = list(structure.get_models())
    if len(models) != 1:
        raise MappingError(f"Expected one PDB model; found {len(models)}")
    chains = list(models[0].get_chains())
    chain_ids = [chain.id for chain in chains]
    if chain_ids != ["A"]:
        raise MappingError(f"Expected only Chain A; found {chain_ids}")

    residues = [residue for residue in chains[0] if is_aa(residue, standard=True)]
    rows: list[dict] = []
    for position, residue in enumerate(residues, start=1):
        aa = seq1(residue.resname, undef_code="X")
        missing = [atom for atom in ("N", "CA", "C") if atom not in residue]
        if missing:
            raise MappingError(
                f"Project position {position} is missing backbone atoms: {missing}"
            )
        if residue.id != (" ", position, " "):
            raise MappingError(
                f"PDB residue ID {residue.id!r} does not equal project position {position}"
            )
        rows.append(
            {
                "position": position,
                "aa": aa,
                "pdb_chain": "A",
                "pdb_resseq": residue.id[1],
            }
        )

    pdb_sequence = "".join(row["aa"] for row in rows)
    exact = len(rows) == EXPECTED_LENGTH and pdb_sequence == reference
    qc = {
        "model_count": len(models),
        "chain_ids": chain_ids,
        "residue_count": len(rows),
        "sequence_exact_match": exact,
        "analysis_allowed": exact,
    }
    if not exact:
        raise MappingError(
            "PDB sequence does not exactly match the frozen 181-residue project reference"
        )
    return rows, qc


def build_official_command(
    *,
    python_executable: str | Path,
    repo_dir: str | Path,
    pdb_path: str | Path,
    out_dir: str | Path,
    models_dir: str | Path,
    cpu_only: bool = False,
) -> list[str]:
    command = [
        str(python_executable),
        str(Path(repo_dir) / "discotope3" / "main.py"),
        "--pdb_or_zip_file",
        str(pdb_path),
        "--struc_type",
        "alphafold",
        "--out_dir",
        str(out_dir),
        "--models_dir",
        str(models_dir),
    ]
    if cpu_only:
        command.append("--cpu_only")
    return command


def read_official_rows(path: str | Path) -> tuple[list[dict], list[str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        missing = [column for column in REQUIRED_OFFICIAL_COLUMNS if column not in fieldnames]
        if missing:
            raise MappingError(f"Official DiscoTope CSV is missing columns: {missing}")
        rows = list(reader)
    return rows, fieldnames


def find_official_csv(
    run_root: str | Path, input_stem: str, chain: str = "A"
) -> Path:
    expected_name = f"{input_stem}_{chain}_discotope3.csv"
    matches = list(Path(run_root).rglob(expected_name))
    if len(matches) != 1:
        raise MappingError(
            f"Expected exactly one official {expected_name}; found {len(matches)}"
        )
    return matches[0]


def preserve_raw_directory(source: str | Path, destination: str | Path) -> Path:
    source = Path(source)
    destination = Path(destination)
    if not source.is_dir():
        raise FileNotFoundError(f"Official raw directory does not exist: {source}")
    if destination.exists():
        raise FileExistsError(f"Raw preservation destination already exists: {destination}")
    shutil.copytree(source, destination)
    return destination


def _finite_float(row: dict, column: str, position: int) -> float:
    try:
        value = float(row[column])
    except (TypeError, ValueError) as error:
        raise MappingError(
            f"Non-numeric {column} at project position {position}: {row[column]!r}"
        ) from error
    if not math.isfinite(value):
        raise MappingError(f"Non-finite {column} at project position {position}")
    return value


def map_official_csv(
    raw_csv: str | Path, reference: str, expected_chain: str = "A"
) -> tuple[list[dict], list[str]]:
    rows, fieldnames = read_official_rows(raw_csv)
    if len(rows) != EXPECTED_LENGTH:
        raise MappingError(
            f"Official DiscoTope output must contain 181 rows; found {len(rows)}"
        )

    mapped: list[dict] = []
    for position, (row, expected_aa) in enumerate(zip(rows, reference), start=1):
        try:
            res_id = int(row["res_id"])
        except ValueError as error:
            raise MappingError(f"Non-integer res_id at output row {position}") from error
        chain = row.get("chain") or row["pdb"].rsplit("_", 1)[-1]
        if chain != expected_chain:
            raise MappingError(
                f"Expected Chain {expected_chain} at project position {position}; found {chain}"
            )
        if res_id != position:
            raise MappingError(
                f"Official res_id {res_id} does not equal project position {position}"
            )
        if row["residue"] != expected_aa:
            raise MappingError(
                f"Residue mismatch at project position {position}: "
                f"expected {expected_aa}, found {row['residue']}"
            )
        if int(row["length"]) != EXPECTED_LENGTH:
            raise MappingError(
                f"Official chain length is not 181 at project position {position}"
            )
        if int(row["alphafold_struc_flag"]) != 1:
            raise MappingError(
                f"AlphaFold structure flag is not 1 at project position {position}"
            )
        mapped.append(
            {
                "position": position,
                "aa": expected_aa,
                "pdb_chain": expected_chain,
                "pdb_resseq": res_id,
                "discotope_score": _finite_float(
                    row, "DiscoTope-3.0_score", position
                ),
                "discotope_calibrated_score": row.get("calibrated_score", ""),
                "discotope_epitope": row.get("epitope", ""),
                "discotope_rsa": _finite_float(row, "rsa", position),
                "discotope_plddt": _finite_float(row, "pLDDTs", position),
                "alphafold_struc_flag": 1,
            }
        )
    return mapped, fieldnames


def _strict_bool(value: object, *, column: str, position: int) -> bool:
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise MappingError(
        f"Non-boolean {column} at project position {position}: {value!r}"
    )


def map_official_csv_phase1b(
    raw_csv: str | Path, reference: str, expected_chain: str = "A"
) -> tuple[list[dict], list[str]]:
    """Map the restored official output to the frozen project schema.

    The uncalibrated official ``DiscoTope-3.0_score`` is retained as
    ``DiscoTope_raw`` for the preregistered percentile transform. The official
    calibrated score and threshold classification are preserved separately.
    """

    official_rows, fieldnames = read_official_rows(raw_csv)
    missing = [
        column
        for column in PHASE1B_ADDITIONAL_OFFICIAL_COLUMNS
        if column not in fieldnames
    ]
    if missing:
        raise MappingError(
            f"Restored official DiscoTope CSV is missing columns: {missing}"
        )

    base_rows, _ = map_official_csv(raw_csv, reference, expected_chain)
    restored: list[dict] = []
    for position, (base, official) in enumerate(
        zip(base_rows, official_rows), start=1
    ):
        restored.append(
            {
                "position": position,
                "AA": base["aa"],
                "DiscoTope_raw": base["discotope_score"],
                "DiscoTope_calibrated": _finite_float(
                    official, "calibrated_score", position
                ),
                "DiscoTope_prediction": _strict_bool(
                    official["epitope"], column="epitope", position=position
                ),
                "DiscoTope_RSA": base["discotope_rsa"],
                "DiscoTope_pLDDT": base["discotope_plddt"],
                "pdb": official["pdb"],
                "chain": expected_chain,
                "pdb_resseq": base["pdb_resseq"],
                "length": int(official["length"]),
                "alphafold_struc_flag": base["alphafold_struc_flag"],
            }
        )
    return restored, fieldnames
