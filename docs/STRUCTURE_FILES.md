# Structure file specification

All coordinates are in Å. In complex PDB files, chain A is IL-24 antigen, H is the antibody heavy chain and L is the antibody light chain. Native uses project residues 27–181; Immunogen uses 27–160. Model-local residue numbering is mapped to project numbering by `data/il24_numbering_map.csv`.

Names identify route, MSA condition and seed. Representative files are unrelaxed Boltz models unless their filename or `results/relaxed_geometry_qc.csv` identifies restrained relaxation. N74-NAG sensitivity structures contain one core GlcNAc only; ordinary structures do not contain that ligand. Relaxed representatives were minimized with restrained backbone heavy atoms; their coordinate RMSD and clash changes are recorded in the QC table. Hydrogen atoms and protonation are those produced by the generating tool and are not experimentally resolved.

Raw Boltz coordinates can contain severe non-bonded clashes. Every PDB in this package is a computational model, not a crystal structure. Candidate interfaces and Phase-1 hotspot overlaps are computational hypotheses, not experimentally confirmed epitopes.
