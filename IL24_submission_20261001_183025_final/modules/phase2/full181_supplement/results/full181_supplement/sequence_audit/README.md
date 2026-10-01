# Full181 sequence and coordinate audit

The uploaded FASTA is exactly 181 amino acids, with line lengths 60/60/60/1,
H180 and L181. Direct sequence comparison agrees with the frozen Full181 FASTA
and every ATOM residue of reference PDB chain A. Native155 is project 27–181;
Immunogen134 is project 27–160. Project 74–76 is NVS. Full local 48 is F and
must not be treated as the old constructs' N74 glycosylation coordinate.

VH119/VL106 agree with both antibody FASTA files after case normalization, all
225 IMGT rows, all 14 CDR/FR ranges, and all 15 available representative PDBs.
The older `data/antibody_numbering.csv` contains inconsistent sequence ranges;
its discrepancies are recorded in `sequence_audit.json` and it is not used.
The frozen `antibody_imgt_numbering.csv` and `cdr_annotation.csv` are consistent.

This audit uses direct sequence equality, explicit residue positions, and path,
size, and modification-time metadata. No digest was computed. The 1–26 boundary
is historical project annotation, not a SignalP or DeepSig measurement. The 15
representative structures do not stand in for 60 original ensemble outputs.

`numbering_map_full181.csv` has 181 rows and includes the old constructs' local
positions, with NA outside a construct. `numbering_map_all_constructs.csv` has
470 rows, preserving separate Full181, Native155, and Immunogen134 maps. These
are input contracts, so their condition is `input_contract` and seed is `NA`.

## Shared helper interface

`src/full181_contract.py` has no import-time writes and requires only Python's
standard library. Supported construct names are Full181, Native155, Immunogen134
and case-insensitive aliases full, native, immunogen.

- `read_constructs(root)` returns validated sequence strings keyed by the three
  canonical names.
- `read_antibodies(root)` returns validated H/L sequence strings.
- `build_mapping(sequence, construct)` validates the entire input sequence and
  returns local_position, project_position, amino_acid rows.
- `build_residue_map(chain, construct, require_complete=True)` accepts a Bio.PDB
  chain or an iterable of `ResidueRecord(id, resname)`. It returns
  `{"residue_map": {residue.id: row}, "metadata": {...}}`. Row keys are
  local_position, project_position, amino_acid, pdb_resseq, and pdb_insertion.
- `build_structure_mapping(pdb_path, construct, root=None, require_complete=True)`
  returns a dictionary keyed A/H/L. Each value has the same residue_map/metadata
  shape. H/L project_position is None. Exact A/H/L protein-chain membership and
  the frozen antigen/antibody sequences are checked.
- `read_pdb_chains(path)` reads ATOM residues in the first PDB model and returns
  per-chain `ResidueRecord` lists, excluding all HETATM records.

Renumbering and insertion-code labels are mapped from residue sequence order,
not arithmetic on PDB residue numbers. Missing residues fail model acceptance
by default. Optional partial mapping accepts only a unique exact subsequence;
ambiguous deletions, substitutions, noncanonical residues, duplicate residue
IDs, and sequence insertions raise ValueError. No deletion position is guessed.

## Verification and reproducibility

Run `python -B tests/test_full181_mapping.py` from the release root. All 9
coordinate tests pass; the exact final output is in `mapping_tests.log`. These
tests check numbering, terminal preservation, sequon positions, changed/truncated
sequence rejection, insertion labels, missing-residue ambiguity, and preserved
legacy antibody/antigen mapping. They do not execute a prediction or smoke test.

Run `python -B src/full181_contract.py --uploaded-fasta <uploaded FASTA path>` to
repeat the audit. Outputs remain under this isolated directory. The exact run
command and exit status are in `sequence_audit_run.log`.

Development verification exposed one incorrect test fixture: the initially
written expected Full local-48 residue was W. Independent rereading of the
uploaded/frozen sequence showed F (project 44–51 is LWEAFWTV), so the fixture was
corrected to F before the final successful run. No biological input was changed.
