#!/usr/bin/env python3
"""Frozen sequence and explicit residue coordinates for the Full181 supplement.

Importing this module performs no I/O. Coordinate maps are keyed by PDB residue
IDs (including insertion codes), never by a presumed PDB-number offset. Missing
residues are accepted only when explicitly requested and uniquely placeable.
"""
import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path

FULL181_SEQUENCE = (
    "MSWGLQILPCLSLILLLWNQVPGLEGQEFRSGSCQVTGVVLPELWEAFWTVKNTVQTQDD"
    "ITSIRLLKPQVLRNVSGAESCYLAHSLLKFYLNTVFKNYHSKIAKFKVLRSFSTLANNFI"
    "VIMSQLQPSKDNSMLPISESAHQRFLLFRRAFKQLDTEVALVKAFGEVDILLTWMQKFYHL"
)
ANTIBODY_SEQUENCES = {
    "H": "EVQLQQSGPELVKPGTSVKVSCKASGYSFTDYNIYWVKQSHGKSLEWIGIDPYNDDTSYNQKFKGKATLTVDKSSSTAFMHLNSLTSEDSAVYYCTRWEITDPWYFDVWGAGTTVTVSS",
    "L": "DIVMTQSHKFMSISVGDRASIMCKASQDVGTAVSWYQQKPGQSPKALTYWASRRHTGVPDRFTGSGSGTDFTLIIGNVQEDLAAYFCQQYSSYPYTFGGGTKLEIK",
}
CONSTRUCTS = {
    "Full181": {"project_start": 1, "project_end": 181, "offset": 0,
                "length": 181, "filename": "IL24_full_1_181.fasta"},
    "Native155": {"project_start": 27, "project_end": 181, "offset": 26,
                  "length": 155, "filename": "IL24_native_27_181.fasta"},
    "Immunogen134": {"project_start": 27, "project_end": 160, "offset": 26,
                     "length": 134, "filename": "IL24_immunogen_27_160.fasta"},
}
_ALIASES = {"full": "Full181", "full181": "Full181", "native": "Native155",
            "native155": "Native155", "immunogen": "Immunogen134",
            "immunogen134": "Immunogen134"}
_AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"))


def construct_name(construct):
    try:
        return _ALIASES[str(construct).lower()]
    except KeyError:
        raise ValueError(f"Unknown construct: {construct}") from None


def expected_sequence(construct):
    spec = CONSTRUCTS[construct_name(construct)]
    return FULL181_SEQUENCE[spec["project_start"]-1:spec["project_end"]]


def sequence_differences(expected, observed):
    return [{"position": i+1,
             "expected": expected[i] if i < len(expected) else None,
             "observed": observed[i] if i < len(observed) else None}
            for i in range(max(len(expected), len(observed)))
            if expected[i:i+1] != observed[i:i+1]]


def _verify_sequence(expected, observed, label):
    differences = sequence_differences(expected, observed)
    if differences:
        raise ValueError(f"{label} sequence mismatch: expected length {len(expected)}, "
                         f"observed {len(observed)}; differences="
                         + json.dumps(differences, separators=(",", ":")))


def read_fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            if not name or name in records:
                raise ValueError(f"Empty or duplicate FASTA header in {path}")
            records[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            records[name] += "".join(line.split()).upper()
    if not records or any(not value for value in records.values()):
        raise ValueError(f"Empty FASTA sequence in {path}")
    return records


def read_constructs(root):
    """Return validated sequences keyed by Full181, Native155, Immunogen134."""
    sequences = {}
    for name, spec in CONSTRUCTS.items():
        records = read_fasta(Path(root) / "data" / spec["filename"])
        if len(records) != 1:
            raise ValueError(f"{name} FASTA must contain exactly one sequence")
        sequence = next(iter(records.values()))
        _verify_sequence(expected_sequence(name), sequence, name)
        sequences[name] = sequence
    return sequences


def read_antibodies(root):
    records = read_fasta(Path(root) / "data/VH_VL_clean.fasta")
    sequences = {}
    for chain, tag in (("H", "VH"), ("L", "VL")):
        matches = [value for key, value in records.items() if tag in key.split()[0]]
        if len(matches) != 1:
            raise ValueError(f"Expected one {tag} FASTA record")
        _verify_sequence(ANTIBODY_SEQUENCES[chain], matches[0], tag)
        sequences[chain] = matches[0]
    if len(records) != 2:
        raise ValueError("Expected precisely VH and VL FASTA records")
    return sequences


def build_mapping(sequence, construct):
    """Validate the entire construct before returning local/project rows."""
    name = construct_name(construct)
    _verify_sequence(expected_sequence(name), sequence, name)
    offset = CONSTRUCTS[name]["offset"]
    return [{"construct": name, "local_position": i,
             "project_position": i + offset, "amino_acid": aa}
            for i, aa in enumerate(sequence, 1)]


@dataclass(frozen=True)
class ResidueRecord:
    """Minimal PDB residue representation, also compatible with Bio.PDB IDs."""
    id: tuple
    resname: str


def read_pdb_chains(path):
    """Read ATOM protein residues from the first PDB model without dependencies."""
    chains = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  "):
            continue
        residue_id = (" ", int(line[22:26]), line[26:27] or " ")
        residue = ResidueRecord(residue_id, line[17:20].strip())
        residues = chains.setdefault(line[21:22], {})
        if residue.id in residues and residues[residue.id].resname != residue.resname:
            raise ValueError(f"Conflicting residue identities at {line[21]}:{residue.id}")
        residues.setdefault(residue.id, residue)
    if not chains:
        raise ValueError(f"No ATOM chains in {path}")
    return {chain: list(residues.values()) for chain, residues in chains.items()}


def _unique_subsequence_positions(observed, expected):
    # Earliest and latest exact subsequence embeddings coincide iff unique.
    earliest, cursor = [], 0
    for aa in observed:
        position = expected.find(aa, cursor)
        if position < 0:
            raise ValueError("Observed sequence is not an exact subsequence; insertion or substitution")
        earliest.append(position + 1)
        cursor = position + 1
    latest, cursor = [], len(expected)
    for aa in reversed(observed):
        position = expected.rfind(aa, 0, cursor)
        latest.append(position + 1)
        cursor = position
    if earliest != list(reversed(latest)):
        raise ValueError("Residue mapping is ambiguous after missing residues; refusing a guessed map")
    return earliest


def _residue_map(chain, expected, offset, label, require_complete):
    residues = [r for r in chain if r.id[0] == " "]
    keys = [tuple(r.id) for r in residues]
    if len(set(keys)) != len(keys):
        raise ValueError(f"Duplicate residue IDs in {label}")
    try:
        observed = "".join(_AA3[r.resname.strip()] for r in residues)
    except KeyError as exc:
        raise ValueError(f"Noncanonical ATOM residue in {label}: {exc}") from None
    if not observed:
        raise ValueError(f"No protein residues in {label}")
    if observed == expected:
        positions = list(range(1, len(expected)+1))
    elif len(observed) < len(expected):
        if require_complete:
            raise ValueError(f"{label} incomplete: {len(observed)}/{len(expected)} residues")
        positions = _unique_subsequence_positions(observed, expected)
    else:
        _verify_sequence(expected, observed, label)
    rows = {key: {"local_position": local,
                  "project_position": local + offset if offset is not None else None,
                  "amino_acid": aa, "pdb_resseq": key[1], "pdb_insertion": key[2].strip()}
            for key, local, aa in zip(keys, positions, observed)}
    return {"residue_map": rows, "metadata": {
        "label": label, "expected_length": len(expected), "observed_length": len(observed),
        "observed_sequence": observed, "complete_sequence_match": observed == expected,
        "missing_local_positions": sorted(set(range(1, len(expected)+1))-set(positions)),
        "has_insertion_codes": any(key[2].strip() for key in keys),
        "local_numbering_contiguous": keys == [(" ", i, " ") for i in range(1, len(expected)+1)],
        "mapping_basis": "exact_sequence" if observed == expected else "unique_exact_subsequence",
    }}


def build_residue_map(chain, construct, require_complete=True):
    """Map a Bio.PDB chain or ResidueRecord iterable to verified project positions.

    Return {'residue_map': {residue.id: row}, 'metadata': {...}}. By default,
    incomplete sequences fail acceptance. To examine a partial structure, set
    require_complete=False; any ambiguous deletion mapping still raises ValueError.
    """
    name = construct_name(construct)
    return _residue_map(chain, expected_sequence(name), CONSTRUCTS[name]["offset"],
                        name, require_complete)


def build_structure_mapping(pdb, construct, root=None, require_complete=True):
    """Return A/H/L residue maps after checking the frozen antigen and antibodies."""
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    read_constructs(root)
    antibody = read_antibodies(root)
    chains = read_pdb_chains(pdb)
    if set(chains) != {"A", "H", "L"}:
        raise ValueError(f"Expected protein chains A/H/L; found {sorted(chains)}")
    result = {"A": build_residue_map(chains["A"], construct, require_complete)}
    for chain in ("H", "L"):
        result[chain] = _residue_map(chains[chain], antibody[chain], None, chain, require_complete)
    return result


def file_metadata(path):
    path = Path(path)
    stat = path.stat()
    return {"source_path": str(path.resolve()), "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "verification_method": "direct sequence/content comparison; filesystem metadata"}


def audit_sequences(root, uploaded_fasta, output):
    """Write isolated audit and numbering tables; never modify legacy inputs."""
    root, uploaded_fasta, output = Path(root), Path(uploaded_fasta), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    sequences, antibody = read_constructs(root), read_antibodies(root)
    upload_records = read_fasta(uploaded_fasta)
    if len(upload_records) != 1:
        raise ValueError("Uploaded FASTA must contain one record")
    upload_sequence = next(iter(upload_records.values()))
    _verify_sequence(FULL181_SEQUENCE, upload_sequence, "uploaded Full181")
    reference_path = root / "data/AF-Q925S4-F1-model_v6.pdb"
    reference = build_residue_map(read_pdb_chains(reference_path)["A"], "Full181")
    imgt_path, cdr_path = root / "data/antibody_imgt_numbering.csv", root / "data/cdr_annotation.csv"
    with imgt_path.open(encoding="utf-8-sig", newline="") as handle:
        imgt = list(csv.DictReader(handle))
    with cdr_path.open(encoding="utf-8-sig", newline="") as handle:
        cdr = list(csv.DictReader(handle))
    annotation_checks = []
    for chain, label in (("H", "VH"), ("L", "VL")):
        records = sorted((r for r in imgt if r["chain"] == label), key=lambda r: int(r["sequence_index"]))
        _verify_sequence(antibody[chain], "".join(r["amino_acid"] for r in records), label + " IMGT")
        if [int(r["sequence_index"]) for r in records] != list(range(1, len(antibody[chain])+1)):
            raise ValueError(f"{label} IMGT sequence indices are not complete and unique")
        coverage = []
        for region in (r for r in cdr if r["chain"] == label):
            start, end = int(region["sequence_start"]), int(region["sequence_end"])
            _verify_sequence(antibody[chain][start-1:end], region["sequence"], label + " " + region["region"])
            if int(region["length"]) != end-start+1:
                raise ValueError(f"{label} {region['region']} length mismatch")
            if any(r["region"] != region["region"] for r in records[start-1:end]):
                raise ValueError(f"{label} IMGT/CDR region inconsistency")
            coverage.extend(range(start, end+1))
        if sorted(coverage) != list(range(1, len(antibody[chain])+1)):
            raise ValueError(f"{label} CDR/FR ranges must cover each residue once")
        annotation_checks.append({"chain": chain, "length": len(antibody[chain]),
                                  "imgt_sequence_and_regions_match": True, "cdr_coverage_complete": True})
    raw_ab = read_fasta(root / "data/VH_VL.fasta")
    for chain, label in (("H", "VH"), ("L", "VL")):
        _verify_sequence(antibody[chain], next(v for k, v in raw_ab.items() if label in k), label + " original FASTA")
    legacy_numbering_path = root / "data/il24_numbering_map.csv"
    with legacy_numbering_path.open(encoding="utf-8-sig", newline="") as handle:
        old_map = list(csv.DictReader(handle))
    if len(old_map) != 181:
        raise ValueError("Legacy IL24 numbering table does not have 181 rows")
    for i, row in enumerate(old_map, 1):
        expected = [str(i), FULL181_SEQUENCE[i-1], str(i-26) if i >= 27 else "",
                    str(i-26) if 27 <= i <= 160 else ""]
        actual = [row[k] for k in ("project_residue_number", "amino_acid", "native_27_181_index", "immunogen_27_160_index")]
        if actual != expected:
            raise ValueError(f"Legacy IL24 numbering mismatch at project position {i}: {actual}")
    legacy_structures = []
    for path in sorted((root / "results").glob("*.pdb")):
        if path.name.startswith("native_") or path.name.startswith("immunogen_"):
            construct = "Native155" if path.name.startswith("native_") else "Immunogen134"
            mapped = build_structure_mapping(path, construct, root)
            legacy_structures.append({**file_metadata(path), "construct": construct,
                                      "all_chain_sequences_match": True,
                                      "chain_metadata": {key: value["metadata"] for key, value in mapped.items()}})
    deprecated_path = root / "data/antibody_numbering.csv"
    deprecated_mismatches = []
    with deprecated_path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            chain = "H" if row["chain"] == "VH" else "L"
            expected = antibody[chain][int(row["start"])-1:int(row["end"])]
            if expected != row["sequence"]:
                deprecated_mismatches.append({"chain": row["chain"], "region": row["region"],
                                              "stored_sequence": row["sequence"], "sequence_at_stored_range": expected})
    sources = [uploaded_fasta, reference_path, imgt_path, cdr_path, legacy_numbering_path,
               root / "data/VH_VL_clean.fasta", root / "data/VH_VL.fasta", deprecated_path]
    sources.extend(root / "data" / spec["filename"] for spec in CONSTRUCTS.values())
    rows = []
    for name, sequence in sequences.items():
        source = str((root / "data" / CONSTRUCTS[name]["filename"]).resolve())
        rows.extend({**row, "condition": "input_contract", "seed": "NA", "source_path": source}
                    for row in build_mapping(sequence, name))
    with (output / "numbering_map_all_constructs.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    full_rows = [{**row,
                  "native_local_position": row["project_position"]-26 if row["project_position"] >= 27 else "NA",
                  "immunogen_local_position": row["project_position"]-26 if 27 <= row["project_position"] <= 160 else "NA"}
                 for row in rows if row["construct"] == "Full181"]
    with (output / "numbering_map_full181.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(full_rows[0]))
        writer.writeheader()
        writer.writerows(full_rows)
    result = {
        "audited_at_utc": datetime.now(timezone.utc).isoformat(), "status": "passed",
        "verification_method": "direct sequence equality and explicit coordinate checks; no digest computation",
        "uploaded_header_metadata_only": list(upload_records)[0],
        "uploaded_sequence_line_lengths": [len(line.strip()) for line in uploaded_fasta.read_text().splitlines()
                                           if line.strip() and not line.startswith(">")],
        "uploaded_length": len(upload_sequence), "uploaded_equals_legacy_full": True,
        "full181_sequence": FULL181_SEQUENCE, "project_180": FULL181_SEQUENCE[179],
        "project_181": FULL181_SEQUENCE[180], "project_74_76": FULL181_SEQUENCE[73:76],
        "native_equals_full_27_181": True, "immunogen_equals_full_27_160": True,
        "reference_chain_A": reference["metadata"], "constructs": CONSTRUCTS,
        "antibody_checks": annotation_checks, "legacy_numbering_map_matches": True,
        "legacy_representative_structures": legacy_structures,
        "legacy_representative_count": len(legacy_structures),
        "deprecated_antibody_numbering_csv": {
            "source_path": str(deprecated_path.resolve()), "used_for_analysis": False,
            "reason": "Stored ranges do not consistently reproduce the frozen antibody sequence",
            "mismatches": deprecated_mismatches},
        "sources": [file_metadata(path) for path in sources],
        "boundary_provenance": "Project positions 1-26 are a historical project definition, not SignalP/DeepSig output",
    }
    (output / "sequence_audit.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--uploaded-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    destination = args.output or args.root / "results/full181_supplement/sequence_audit"
    result = audit_sequences(args.root, args.uploaded_fasta, destination)
    print(json.dumps({"status": result["status"], "length": result["uploaded_length"],
                      "legacy_representative_count": result["legacy_representative_count"],
                      "output": str(destination)}))
