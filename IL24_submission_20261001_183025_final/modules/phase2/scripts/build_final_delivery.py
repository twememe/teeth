#!/usr/bin/env python3
"""Build the auditable Phase 2 result archive and SHA-256 manifests."""

from __future__ import annotations

import csv
import hashlib
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DELIVERY = RESULTS / "17_delivery"
MANIFEST = DELIVERY / "final_artifact_manifest_sha256.csv"
ARCHIVE = RESULTS / "phase2_final_results.tar.gz"
SIDECAR = RESULTS / "phase2_final_results.sha256"


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def selected_files() -> list[Path]:
    roots = [
        ROOT / "data" / "IL24_full_1_181.fasta",
        ROOT / "data" / "IL24_native_27_181.fasta",
        ROOT / "data" / "IL24_immunogen_27_160.fasta",
        ROOT / "data" / "VH_VL_clean.fasta",
        ROOT / "data" / "il24_numbering_map.csv",
        ROOT / "data" / "phase2_input_manifest_sha256.csv",
        ROOT / "data" / "antibody_imgt_numbering.csv",
        ROOT / "data" / "cdr_annotation.csv",
        ROOT / "data" / "msa",
        ROOT / "scripts",
        RESULTS / "11_antibody_structure",
        RESULTS / "12_complex_prediction" / "formal",
        RESULTS / "12_complex_prediction" / "formal_msa",
        RESULTS / "12_complex_prediction" / "msa_sensitivity",
        RESULTS / "12_complex_prediction" / "glyco_sensitivity",
        RESULTS / "12_complex_prediction" / "glyco_schema_validation",
        RESULTS / "12_complex_prediction" / "boltz_inputs",
        RESULTS / "12_complex_prediction" / "phase1_prior_deduplicated.csv",
        RESULTS / "12_complex_prediction" / "sequences_info.json",
        RESULTS / "12_complex_prediction" / "step12_input_manifest_sha256.csv",
        RESULTS / "13_interface_analysis",
        RESULTS / "14_quality",
        RESULTS / "15_affinity",
        RESULTS / "16_mechanism",
        DELIVERY,
        ROOT / "logs",
    ]
    files: set[Path] = set()
    for item in roots:
        if item.is_file():
            files.add(item)
        elif item.is_dir():
            files.update(
                path
                for path in item.rglob("*")
                if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc"
            )
    files.discard(MANIFEST)
    files.discard(ARCHIVE)
    files.discard(SIDECAR)
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def main() -> None:
    DELIVERY.mkdir(parents=True, exist_ok=True)
    files = selected_files()
    rows = [
        {
            "file": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
        for path in files
    ]
    for relative, size, sha256 in (
        ("models/boltz/boltz1_conf.ckpt", 3595352714, "fea245d912c570ec117b2277c2719f312a6fc109c07b6f6ef741690ee775c2f5"),
        ("models/boltz/ccd.pkl", 345859128, "2d3b2f03a3c5665944adba51e33263511e51b21c9cd05d902f9c4b7c1e58d2f4"),
    ):
        rows.append({"file": relative, "bytes": size, "sha256": sha256})
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("file", "bytes", "sha256"))
        writer.writeheader()
        writer.writerows(rows)

    archive_files = files + [MANIFEST]
    with tarfile.open(ARCHIVE, "w:gz", compresslevel=6) as archive:
        for path in archive_files:
            archive.add(path, arcname=path.relative_to(ROOT).as_posix(), recursive=False)
    archive_sha = digest(ARCHIVE)
    SIDECAR.write_bytes(f"{archive_sha}  {ARCHIVE.name}\n".encode("ascii"))
    print(f"Manifest entries: {len(rows)}")
    print(f"Archive: {ARCHIVE} ({ARCHIVE.stat().st_size} bytes)")
    print(f"SHA256: {archive_sha}")


if __name__ == "__main__":
    main()
