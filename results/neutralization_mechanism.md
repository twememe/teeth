# IL-24 neutralization mechanism assessment

## Evidence available

- Native 27-181 models show a >=50% consensus interface around project residues 57-65, 98-121 and 163.
- Immunogen 27-160 models show a different >=50% consensus interface around residue 62 and 133-153.
- The Native/Immunogen consensus-interface Jaccard index is 0.031.
- Most modeled antibody contacts map to CDR residues, but complex ipTM is only moderate (Native mean 0.505; Immunogen mean 0.535).
- In a matched three-seed, 50-step sensitivity test, ColabFold MSAs increased Native/Immunogen consensus-interface Jaccard from 0.000 to 0.632 and increased complex pLDDT, showing that the empty-MSA formal interface is not robust to MSA choice.
- In the completed 15-seed/200-step formal MSA ensembles, Native/Immunogen consensus-interface Jaccard is 0.667 versus 0.031 in the empty-MSA ensembles. The shared MSA consensus centers on project residues 135-150.
- Strict shared 27-160 backbone alignment gives a median Native/Immunogen antibody-pose C-alpha RMSD of 42.24 A for the formal empty-MSA ensembles and 16.91 A for the MSA sensitivity ensembles; the MSA result remains heterogeneous (8.33-63.15 A across all seed pairs).
- N74 is the only N-X-S/T sequon. It contacts the antibody directly in 2/3 MSA Native models but 0/30 formal empty-MSA models, so unmodeled glycosylation can materially alter the inferred interface.
- A matched single-NAG model caused seed-dependent interface rearrangement and a 1.51 A NAG-antibody clash in Native seed 2, excluding that pose as sterically plausible without further relaxation/reorientation.
- Five shortlisted representatives were restrained-minimized; all `<1.8 A` severe clashes were removed with only 0.108-0.191 A backbone RMSD.

## Mechanistic conclusion

No reliable IL-24-receptor complex or experimentally established receptor-contact table is present inside the permitted Phase 2 workspace. Direct receptor competition, steric obstruction and allosteric interference therefore cannot be distinguished computationally from the current data.

The formal MSA consensus at 135-150 supports the broader Phase 1B 126-157 region as an experimental candidate, not as a proven epitope. Phase 1B is not ground truth. The disagreement with empty-MSA models, remaining pose heterogeneity and N74 glycosylation sensitivity prevent a definitive binding-mode claim even though representative-model clashes can be repaired by restrained minimization.

## Recommended validation

1. Mutagenesis or peptide/alanine scanning across 57-65, 117-121 and 133-153.
2. Competition binding against IL-20RA/IL-20RB and IL-22RA1/IL-20RB receptor complexes.
3. SPR/BLI affinity measurements for native mature IL-24 and the 27-160 immunogen.
4. Model glycosylated N74 explicitly and compare with unglycosylated controls.
5. Expand the MSA condition beyond three seeds, relax shortlisted interfaces, and cross-check with a second prediction engine.
