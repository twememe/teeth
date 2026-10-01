# Data and provenance

Project sequences and antibody inputs are stored in this directory. External structure, papers and ELISA data are in `external_materials/` with SHA-256 checksums. RCSB 6DF3 is public structural data from https://www.rcsb.org/structure/6DF3; the two papers are retained as supplied reference documents; ELISA CSV is project-generated experimental data.

The input manifest records hashes for frozen inputs. The workflow does not use a hidden evaluation set. Native and Immunogen are analyzed as separate constructs; MSA files are retained with a manifest. Results are hypothesis-generation outputs and require experimental validation.
