#!/usr/bin/env python3
"""Generate Boltz YAML inputs from the frozen Phase 2 FASTA files."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "12_complex_prediction" / "boltz_inputs"


def read_fasta(path):
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name:
            records[name] += line.strip().upper()
    return records


ab = read_fasta(DATA / "VH_VL_clean.fasta")
vh = next(v for k, v in ab.items() if "VH" in k)
vl = next(v for k, v in ab.items() if "VL" in k)


def antigen(filename):
    return next(iter(read_fasta(DATA / filename).values()))


def yaml_text(seq, msa_paths=None):
    msa_paths = msa_paths or {"A": "empty", "H": "empty", "L": "empty"}
    return f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}
      msa: {msa_paths['A']}
  - protein:
      id: H
      sequence: {vh}
      msa: {msa_paths['H']}
  - protein:
      id: L
      sequence: {vl}
      msa: {msa_paths['L']}
"""


def glyco_yaml_text(seq, target):
    msa_paths = {chain: f"data/msa/{target}/{chain}.csv" for chain in "AHL"}
    return yaml_text(seq, msa_paths) + """  - ligand:
      id: G
      ccd: NAG
constraints:
  - bond:
      atom1: [A, 48, ND2]
      atom2: [G, 1, C1]
"""


OUT.mkdir(parents=True, exist_ok=True)
(OUT / "native_blind.yaml").write_text(yaml_text(antigen("IL24_native_27_181.fasta")))
(OUT / "native_phase1_informed.yaml").write_text(yaml_text(antigen("IL24_native_27_181.fasta")))
(OUT / "immunogen_blind.yaml").write_text(yaml_text(antigen("IL24_immunogen_27_160.fasta")))
(OUT / "native_blind_msa.yaml").write_text(yaml_text(
    antigen("IL24_native_27_181.fasta"),
    {chain: f"data/msa/native/{chain}.csv" for chain in "AHL"},
))
(OUT / "immunogen_blind_msa.yaml").write_text(yaml_text(
    antigen("IL24_immunogen_27_160.fasta"),
    {chain: f"data/msa/immunogen/{chain}.csv" for chain in "AHL"},
))
(OUT / "native_blind_msa_n74_nag.yaml").write_text(glyco_yaml_text(
    antigen("IL24_native_27_181.fasta"), "native",
))
(OUT / "immunogen_blind_msa_n74_nag.yaml").write_text(glyco_yaml_text(
    antigen("IL24_immunogen_27_160.fasta"), "immunogen",
))
print(f"Wrote 7 inputs to {OUT}; VH={len(vh)}, VL={len(vl)}")
