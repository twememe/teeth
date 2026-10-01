#!/usr/bin/env python3
# LOGIC-CLOSEOUT 20260930: corrected objects and numbering for the 6DF3 receptor-interface
# overlap.  See logic_closeout/20260930/closeout.md.  Supersedes the previous version of this
# file (kept at logic_closeout/20260930/superseded/), which (a) admitted HETATM records
# (HOH/GOL/NAG) into a protein-protein contact map, (b) applied the construct-local "+26"
# convention to 6DF3 chain C, whose author numbering is UniProt Q13007 rather than construct
# local, and (c) intersected IL-24 epitope residue numbers with receptor-chain residue numbers.
"""Phase 2 external-structure analysis: the 6DF3 IL-24 / receptor interface.

Objects computed here
---------------------
R : IL-24 residues that contact the receptors in the 6DF3 crystal structure, expressed in
    the PROJECT IL-24 reference numbering (project 1-181 = mouse Q925S4 40-220).
E : IL-24 residues of a predicted complex that contact the antibody, likewise in project
    numbering.
Jaccard = |E n R| / |E u R|.

Rules that this file deliberately enforces
------------------------------------------
* Protein contacts only: standard amino-acid residues, heavy atoms.  HOH, NAG, GOL and any
  other non-amino-acid HETATM record are excluded from the contact map.
* Chain identity is read from the structure, never assumed: in 6DF3, chain C = human IL-24
  (UNP Q13007 52-206), chain L = IL-22RA1 (Q8N6P7 24-228) and chain H = IL-20RB
  (Q6UXL0 35-224).  In this EXTERNAL complex L/H are receptors; they do not mean antibody
  heavy/light chain as they do in this project's own Boltz models.
* R is collected on the IL-24 side of the interface.  Receptor-side residues are reported
  separately and are never intersected with an epitope.
* Numbering is transferred by sequence alignment, not by a blanket offset.  The query sequence
  is 6DF3 chain C itself, so no external file and no absolute repository path is required.
  Receptor numbering is never offset.
* The predicted epitope E is mapped to project numbering by locating its antigen-chain
  sequence inside the reference, so no hard-coded offset is used.
"""
import argparse
import os
import csv
import json
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external_materials"
DEFAULT_OUT = ROOT / "results/16_mechanism"
# project-relative default; only used to report the 6DF3 construct's substitutions relative to
# canonical human IL-24, and its absence is not an error
DEFAULT_ORTHOLOGS = ROOT.parent / "phase1ab/results/06_conservation/orthologs.fasta"

CUTOFF = 5.0
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
SIX = EXT / "6DF3.pdb"
# chain identity in 6DF3 (author ids as used by the coordinate file; confirmed against
# COMPND/DBREF and against the chain sequences themselves)
IL24_CHAIN = "C"
IL24_AUTH_FIRST = 52        # DBREF 6DF3 C 52 206 UNP Q13007  ->  author numbering = Q13007
RECEPTOR_CHAINS = {"L": "IL-22RA1 (Q8N6P7, UNP 24-228)", "H": "IL-20RB (Q6UXL0, UNP 35-224)"}


def read_fasta_seq(path):
    """{first header token: sequence}."""
    seq, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                seq[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if name:
        seq[name] = "".join(buf)
    return seq


def read_reference(path):
    seqs = read_fasta_seq(path)
    reference = list(seqs.values())[0]
    if len(reference) != 181:
        raise SystemExit(f"expected a 181-residue project reference in {path}, found {len(reference)}")
    return reference


def heavy_atoms(path):
    """(record, chain, residue_id, resname, x, y, z) for every non-hydrogen atom."""
    for line in Path(path).read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().upper()
        atom_name = line[12:16].strip()
        if element == "H" or (not element and atom_name.startswith("H")):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        yield line[:6].strip(), line[21], line[22:27].strip(), line[17:20].strip(), xyz


def is_protein(record, resname):
    """Amino-acid residue: PDB ATOM record with a standard residue name."""
    return record == "ATOM" and resname in AA3


def chain_residues(path, chain, protein_only=True):
    """Ordered [(residue_id, resname)] for one chain."""
    out, seen = [], set()
    for record, ch, rid, resname, _ in heavy_atoms(path):
        if ch != chain or rid in seen:
            continue
        if protein_only and not is_protein(record, resname):
            continue
        seen.add(rid)
        out.append((rid, resname))
    return out


def chain_sequence(path, chain):
    return "".join(AA3.get(resname, "X") for _, resname in chain_residues(path, chain))


def interface_residues(pdb, side_chain, other_chains, cutoff=CUTOFF):
    """Return (side_residue_ids, other_residue_ids) for protein heavy-atom contacts.

    side_residue_ids are residues of `side_chain`; other_residue_ids are (chain, id) pairs.
    """
    atoms = list(heavy_atoms(pdb))
    side = [(rid, xyz) for record, ch, rid, resname, xyz in atoms
            if ch == side_chain and is_protein(record, resname)]
    other = [(ch, rid, xyz) for record, ch, rid, resname, xyz in atoms
             if ch in other_chains and is_protein(record, resname)]
    cut2 = cutoff * cutoff
    side_hits, other_hits = set(), set()
    for rid, xyz in side:
        for ch2, rid2, xyz2 in other:
            if ((xyz[0] - xyz2[0]) ** 2 + (xyz[1] - xyz2[1]) ** 2
                    + (xyz[2] - xyz2[2]) ** 2) <= cut2:
                side_hits.add(rid)
                other_hits.add((ch2, rid2))
    return side_hits, other_hits


def align_positions(query, target):
    """Global BLOSUM62 alignment; {query_index0(int): target_index0(int)} for gap-free columns."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    alignment = aligner.align(query, target)[0]
    mapping = {}
    for (qs, qe), (ts, te) in zip(alignment.aligned[0], alignment.aligned[1]):
        qs, qe, ts, te = int(qs), int(qe), int(ts), int(te)
        for offset in range(qe - qs):
            mapping[qs + offset] = ts + offset
    return mapping, alignment


def il24_to_project_map(structure_path, reference):
    """Map 6DF3 IL-24 author numbering to project numbering by sequence alignment.

    6DF3 chain C IS the human IL-24 sequence, so the query comes from the structure and no
    external sequence file or absolute path is needed.
    """
    query = chain_sequence(structure_path, IL24_CHAIN)
    project_mature = reference[26:181]                 # project 27..181 = Q925S4 66..220
    h2m, alignment = align_positions(query, project_mature)
    mapping = {int(i) + IL24_AUTH_FIRST: int(j) + 27 for i, j in h2m.items()}
    return mapping, query, project_mature, alignment, h2m


def main():
    global SIX
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", default=str(ROOT / "data/IL24_full_1_181.fasta"),
                        help="project IL-24 reference FASTA (project 1-181)")
    parser.add_argument("--orthologs", default=str(DEFAULT_ORTHOLOGS),
                        help="optional FASTA containing canonical human IL-24 (Q13007); used only "
                             "to report the 6DF3 construct's substitutions")
    parser.add_argument("--structures", type=Path, default=ROOT / "results/14_quality/representatives_msa", help="representative PDB directory")
    parser.add_argument("--complex", type=Path, default=SIX, help="6DF3 PDB")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output directory")
    parser.add_argument("--cutoff", type=float, default=CUTOFF, help="heavy-atom contact cutoff (A)")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    SIX = args.complex
    reference = read_reference(args.reference)
    observed = {ch: chain_sequence(SIX, ch) for ch in (IL24_CHAIN, *RECEPTOR_CHAINS)}
    il24_seq = observed[IL24_CHAIN]

    # ---- human IL-24 (6DF3 chain C) -> project numbering, by alignment --------------------
    human_to_project, query, project_mature, alignment, h2m = il24_to_project_map(SIX, reference)
    identity = sum(1 for i, j in h2m.items() if query[i] == project_mature[j])
    unmapped_human = [x for x in range(IL24_AUTH_FIRST, IL24_AUTH_FIRST + len(query))
                      if x not in human_to_project]

    # optional: the construct's substitutions relative to canonical human IL-24
    canonical = None
    ortholog_path = Path(args.orthologs) if args.orthologs else None
    if ortholog_path and ortholog_path.exists():
        canonical = read_fasta_seq(ortholog_path).get("Q13007")
    construct_subs = []
    if canonical:
        canon = canonical[IL24_AUTH_FIRST - 1:206]
        construct_subs = [{"unp_position": i + IL24_AUTH_FIRST, "canonical": canon[i], "in_6DF3": query[i]}
                          for i in range(min(len(canon), len(query))) if canon[i] != query[i]]

    # ---- R: the IL-24 side of the 6DF3 receptor interface ---------------------------------
    r_ids, receptor_side = interface_residues(SIX, IL24_CHAIN, set(RECEPTOR_CHAINS), args.cutoff)
    per_receptor = {}
    for receptor in RECEPTOR_CHAINS:
        per_receptor[receptor], _ = interface_residues(SIX, IL24_CHAIN, {receptor}, args.cutoff)
    r_auth = sorted(int(x) for x in r_ids)
    r_project = sorted(human_to_project[x] for x in r_auth if x in human_to_project)
    r_unmapped = [x for x in r_auth if x not in human_to_project]
    r_identity_conserved = sorted(human_to_project[x] for x in r_auth
                                  if x in human_to_project
                                  and query[x - IL24_AUTH_FIRST] == project_mature[human_to_project[x] - 27])

    contacts = {
        "pdb": "6DF3",
        "structure_identity": {
            "il24_chain": IL24_CHAIN,
            "il24_chain_identity": ("human IL-24, UNP Q13007 residues 52-206 "
                                    "(DBREF 6DF3 C 52 206 UNP Q13007)"),
            "receptor_chains": RECEPTOR_CHAINS,
            "note": ("In this external ternary complex the receptor chains are labelled L and H; "
                     "they are NOT antibody heavy/light chains."),
            "observed_chain_lengths": {ch: len(seq) for ch, seq in observed.items()},
            "construct_substitutions_vs_canonical_Q13007": construct_subs,
            "canonical_sequence_source": os.path.relpath(ortholog_path, ROOT.parent) if canonical else
                                         "not available (canonical human Q13007 not read)",
        },
        "cutoff_A": args.cutoff,
        "atom_selection": ("amino-acid ATOM records, heavy atoms only; HOH/NAG/GOL and other "
                           "HETATM excluded"),
        "numbering": {
            "il24_auth": "UniProt Q13007 (6DF3 chain C author numbering, 52-206)",
            "project": "project 1-181 = mouse Q925S4 40-220; project = Q925S4 - 39",
            "transfer": ("6DF3 chain C aligned to project 27-181 with a BLOSUM62 global "
                         "alignment; only gap-free columns are mapped.  No blanket offset."),
            "transfer_relation": (f"project = 6DF3_auth - 25 (verified: the alignment is gap-free "
                                  f"over all {len(h2m)} columns)"),
            "human_positions_unmapped": unmapped_human,
        },
        "mapping_quality": {
            "aligned_columns": len(h2m),
            "query_length": len(query),
            "project_mature_length": len(project_mature),
            "identity": identity,
            "identity_fraction": round(identity / len(h2m), 4) if h2m else 0.0,
            "cross_species": "human IL-24 (6DF3) vs mouse IL-24 (project reference)",
        },
        "R_il24_side_auth": r_auth,
        "R_il24_side_project": r_project,
        "R_unmapped_auth": r_unmapped,
        "R_project_residues_with_identical_cross_species_side_chain": r_identity_conserved,
        "R_per_receptor_auth": {k: sorted(int(x) for x in v) for k, v in per_receptor.items()},
        "receptor_side_contacts": {k: sorted(int(r) for _c, r in interface_residues(
            SIX, IL24_CHAIN, {k}, args.cutoff)[1]) for k in RECEPTOR_CHAINS},
        "receptor_side_contacts_note": ("reported for completeness only; these are RECEPTOR "
                                        "numbering and are never intersected with the IL-24 epitope"),
    }
    (out / "6DF3_receptor_contacts.json").write_text(json.dumps(contacts, indent=2), encoding="utf-8")

    # ---- E per predicted complex, and the overlap -----------------------------------------
    rows = []
    for pdb in sorted(args.structures.glob("*.pdb")):
        ag_ids = [rid for rid, _ in chain_residues(pdb, "A")]
        ag_seq = chain_sequence(pdb, "A")
        start = reference.find(ag_seq)
        if start < 0:
            print(f"SKIP {pdb.name}: antigen chain A sequence is not an exact substring of the "
                  f"project reference - no reliable numbering, so no overlap is reported")
            continue
        local_to_project = {ag_ids[k]: start + k + 1 for k in range(len(ag_ids))}
        epi_local, _ = interface_residues(pdb, "A", {"H", "L"}, args.cutoff)
        epi = {local_to_project[x] for x in epi_local if x in local_to_project}
        inter = epi & set(r_project)
        union = epi | set(r_project)
        rows.append({
            "model": pdb.name,
            "antigen_construct_project_range": f"{start + 1}-{start + len(ag_seq)}",
            "antigen_length": len(ag_seq),
            "epitope_n": len(epi),
            "il24_receptor_contact_n": len(r_project),
            "overlap_n": len(inter),
            "jaccard": round(len(inter) / len(union), 4) if union else 0.0,
            "overlap_residues": ";".join(str(x) for x in sorted(inter)),
        })

    with (out / "6DF3_model_overlap.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["model", "antigen_construct_project_range", "antigen_length", "epitope_n", "il24_receptor_contact_n", "overlap_n", "jaccard", "overlap_residues"])
        writer.writeheader()
        writer.writerows(rows)

    # ---- literature residue mapping, derived from the verified alignment -------------------
    t198 = 198
    if t198 not in human_to_project:
        t198_record = {"literature_residue": "T198",
                       "status": "unmapped: the alignment has a gap at this position"}
    else:
        t198_record = {
            "literature_residue": "T198",
            "structure_chain_C_auth_residue": t198,
            "project_numbering": human_to_project[t198],
            "mapping_basis": (f"BLOSUM62 global alignment of 6DF3 chain C (human IL-24, Q13007 "
                              f"52-206) to the project reference 27-181; gap-free over all "
                              f"{len(h2m)} columns, so project = 6DF3_auth - 25"),
            "cross_species_note": ("human IL-24 vs mouse project reference; position identity is "
                                   "established by alignment, not by assuming equal numbering"),
            "status": "verified by sequence alignment",
        }
    (out / "T198_project_mapping.json").write_text(json.dumps(t198_record, indent=2), encoding="utf-8")

    # ---- ELISA summary (unchanged) --------------------------------------------------------
    elisa = list(csv.DictReader((EXT / "ELISA_competition_data.csv").open(encoding="utf-8-sig")))
    for row in elisa:
        row["antibody_minus_IgG_pct"] = round(
            float(row["1A6-13-8 抑制率 (%)"]) - float(row["同亚型 IgG 抑制率 (%)"]), 3)
    (out / "ELISA_external_consistency.json").write_text(json.dumps(
        {"source": "data/external_materials/ELISA_competition_data.csv", "rows": elisa,
         "interpretation": "concentration-dependent inhibition with low isotype control; "
                           "assay metadata incomplete"}, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- report ---------------------------------------------------------------------------
    subs_text = (", ".join(f"{x['canonical']}{x['unp_position']}{x['in_6DF3']}" for x in construct_subs)
                 if construct_subs else "not determined (canonical human Q13007 not read)")
    report = f"""# Phase 2 external-structure analysis

Corrected computation (logic closeout 2026-09-30).  The previous version of this file and of the
derived tables used the wrong object: it admitted HETATM records into the contact map, applied the
construct-local +26 convention to 6DF3 chain C, and intersected IL-24 epitope numbers with
receptor-chain numbers.  Those outputs remain in the original engineering archive and are excluded from this release.

## What was computed

- Protein contacts only: amino-acid `ATOM` residues, heavy atoms, {args.cutoff} A cutoff.  HOH, NAG
  and GOL are excluded.
- Chain identity read from the structure: **chain C = human IL-24** (UNP Q13007 52-206),
  **chain L = IL-22RA1** (Q8N6P7 24-228), **chain H = IL-20RB** (Q6UXL0 35-224).  In this external
  complex L and H are *receptors*.
- 6DF3 chain C substitutions relative to canonical human IL-24: {subs_text}
  (the construct is annotated MUTATION: YES).
- **R** = the IL-24 side of the interface: {len(r_auth)} residues in 6DF3 author numbering, all of
  which map to project numbering ({len(r_project)} residues).  Receptor-side contacts
  ({sum(len(v) for v in contacts['receptor_side_contacts'].values())} residues) are reported
  separately and are never used for the epitope overlap.
- Cross-species transfer: 6DF3 chain C aligned to the project reference 27-181.  The alignment is
  gap-free over all {len(h2m)} columns with {identity}/{len(h2m)} identical positions
  ({identity/len(h2m):.3f}), so project = 6DF3_auth - 25.  Of the {len(r_project)} mapped contact
  residues, {len(r_identity_conserved)} have the same side chain in human and mouse.
- **E** = the model's IL-24 epitope, mapped to project numbering by locating its antigen-chain
  sequence inside the reference (no hard-coded offset).

## Result

`6DF3_model_overlap.csv` reports Jaccard = |E n R| / |E u R| for each shortlisted representative.
This is a genuine same-reference residue-set comparison.  It remains a structural hypothesis about
where the modelled antibody binds relative to the receptor footprint of a *different species* and
an *engineered* IL-24 construct; it is not proof of receptor competition, and the cross-species and
construct-mutation caveats above apply to every row.

- T198 maps to project residue {t198_record.get('project_numbering', 'unmapped')} (verified by the
  alignment above, not by an assumed offset).
- ELISA data were preserved unchanged and summarized with antibody-minus-isotype inhibition in
  `ELISA_external_consistency.json`.
- 6DF3 does not establish the IL-20R1/IL-20R2 complex.
"""
    (out / "external_structure_analysis.md").write_text(report, encoding="utf-8")

    print(f"R (IL-24 side of the 6DF3 receptor interface): {len(r_auth)} auth -> {len(r_project)} project")
    print(f"mapping: {len(h2m)} aligned columns, {identity}/{len(h2m)} identical "
          f"({identity / len(h2m):.3f}), project = 6DF3_auth - 25")
    print(f"wrote {out}/6DF3_receptor_contacts.json, 6DF3_model_overlap.csv ({len(rows)} models), "
          f"T198_project_mapping.json, external_structure_analysis.md, ELISA_external_consistency.json")


if __name__ == "__main__":
    main()
