# Structural sensitivity and geometry audit

## MSA sensitivity

The controlled comparison used identical Boltz-1 settings (3 recycling steps, 50 sampling steps, one diffusion sample) and seeds 1-3. Only the MSA condition changed. Raw MSA depths were 1,475-11,310 sequences depending on target/chain; Boltz used its default effective cap of 8,192 sequences per chain.

- Native: mean ipTM increased from 0.501 to 0.529, mean complex pLDDT from 0.739 to 0.848, and pairwise epitope Jaccard from 0.359 to 0.590.
- Immunogen: mean ipTM increased from 0.554 to 0.565 and mean complex pLDDT from 0.744 to 0.823; pairwise epitope Jaccard changed only from 0.234 to 0.250.
- Native hotspot coverage increased from 0.188 to 0.247, while Immunogen hotspot coverage decreased from 0.301 to 0.167.
- Native/Immunogen >=50% consensus-interface Jaccard in the shared 27-160 region increased from 0.000 to 0.632.

The MSA condition materially changes the inferred binding interface. Three seeds establish sensitivity but do not replace the 15-seed formal ensembles.

The subsequent full 15-seed/200-step MSA ensembles confirmed the main sensitivity result:

- Native mean ipTM increased from 0.505 to 0.549 and complex pLDDT from 0.743 to 0.851.
- Immunogen mean ipTM increased from 0.535 to 0.544 and complex pLDDT from 0.741 to 0.817.
- Native/Immunogen >=50% consensus-interface Jaccard increased from 0.031 to 0.667; the shared MSA consensus is concentrated at project residues 135, 136, 139, 140, 142, 143, 146 and 150.
- Mean within-route pairwise epitope Jaccard changed from 0.450 to 0.399 for Native and from 0.235 to 0.314 for Immunogen. The MSA ensemble is therefore still multimodal despite stronger cross-construct consensus.

## Native/Immunogen pose

Each Immunogen model was aligned to each Native model using the 402 shared antigen backbone atoms (N, CA and C of local residues 1-134/project residues 27-160). Antibody pose RMSD was then calculated across the 225 matched H/L CA atoms without independently realigning the antibody.

- Formal empty-MSA ensembles (225 cross-pairs): median shared-antigen RMSD 12.30 A; median antibody-pose RMSD 42.24 A.
- MSA sensitivity ensembles (9 cross-pairs): median shared-antigen RMSD 1.97 A; median antibody-pose RMSD 16.91 A, range 8.33-63.15 A.
- Formal MSA ensembles (225 cross-pairs): median shared-antigen RMSD 3.50 A; median antibody-pose RMSD 30.66 A, range 1.96-64.91 A.

MSA improves antigen structural agreement and the median pose comparison, but antibody binding poses remain heterogeneous.

## N-glycosylation sensitivity

N74 (NVS; local residue 48) is the only N-X-S/T sequon in IL-24. A matched three-seed/200-step sensitivity run requested one NAG on ASN ND2 (exploratory): the bond is present in the output topology, but the realised ND2-C1 distance is 1.58-2.20 A against about 1.44-1.46 A for the same bond in the 6DF3 reference, and no constraint-residual check was applied. This is not evidence that the complex is stable in the presence of a glycan, and the N74 distances in results/13_interface_analysis/glycosylation_interface_analysis.csv are distances from the protein N74 residue, not from glycan atoms. Relative to the same unmodified MSA seeds, Native mean ipTM changed by -0.018, complex pLDDT by -0.024 and hotspot coverage by -0.038. Mean paired epitope Jaccard was 0.549 for Native and 0.503 for Immunogen, with strong seed dependence. Native seed 2 placed NAG only 1.51 A from the antibody, a severe steric clash. This demonstrates glycosylation sensitivity but represents only the core GlcNAc, not a full heterogeneous human N-glycan.

## Geometry

No empty-MSA formal model has a CA-defined chain break or peptide C-N outlier. All 66 audited unglycosylated, unrelaxed Boltz structures contain at least one severe nonbonded heavy-atom clash below 1.8 A; the 50-step MSA Immunogen seed-2 model also has one peptide C-N outlier. These coordinates are suitable for ensemble-level hypothesis generation and relative ranking, not atomistic interaction claims. Any shortlisted structure should be repaired/relaxed before detailed energetic or mechanism interpretation.

Five representatives were subsequently minimized with Amber14/GBn2 and 1000 kJ mol-1 nm-2 backbone restraints. All severe interchain and nonlocal intrachain clashes were removed, while backbone coordinate RMSD remained 0.108-0.191 A. Relaxed coordinates are in `results/14_quality/relaxed/`.
