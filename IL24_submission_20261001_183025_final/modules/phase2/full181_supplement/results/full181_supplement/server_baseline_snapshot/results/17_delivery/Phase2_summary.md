# Phase 2 computational summary

Formal models analyzed: 30 empty-MSA + 30 ColabFold-MSA; plus 6 N74-NAG sensitivity models.

## Empty-MSA formal route summaries

### immunogen_blind

- n=15; mean ipTM=0.535 (SD 0.008)
- mean complex pLDDT=0.741
- mean Phase-1 hotspot coverage=0.294
- >=50% epitope consensus residues (project numbering): 62;133;134;135;136;137;139;142;145;146;149;153

### native_blind

- n=15; mean ipTM=0.505 (SD 0.007)
- mean complex pLDDT=0.743
- mean Phase-1 hotspot coverage=0.157
- >=50% epitope consensus residues (project numbering): 57;58;59;62;63;65;98;99;101;103;104;106;107;108;110;111;114;117;118;121;163

## Native versus immunogen

- Jaccard index of >=50% epitope consensus sets: 0.031
- Comparison uses project numbering and the shared 27-160 region; residue 163 in Native is outside the immunogen construct.

## Formal 15-seed MSA ensembles

- native_blind: mean ipTM 0.505 -> 0.549; complex pLDDT 0.743 -> 0.851; mean pairwise epitope Jaccard 0.450 -> 0.399; hotspot coverage 0.157 -> 0.189.
- immunogen_blind: mean ipTM 0.535 -> 0.544; complex pLDDT 0.741 -> 0.817; mean pairwise epitope Jaccard 0.235 -> 0.314; hotspot coverage 0.294 -> 0.170.
- Native/Immunogen >=50% consensus Jaccard in the shared region is 0.667 with formal MSA, versus 0.031 with empty MSA.
- The formal-MSA consensus is concentrated at project residues 132-150, but seed-level pose heterogeneity remains substantial.

## Matched MSA sensitivity analysis

Three matched seeds per route used 50 sampling steps for both empty-MSA and ColabFold-MSA conditions.

- native_blind: ipTM 0.501 -> 0.529; complex pLDDT 0.739 -> 0.848; pairwise epitope Jaccard 0.359 -> 0.590; hotspot coverage 0.188 -> 0.247.
- immunogen_blind: ipTM 0.554 -> 0.565; complex pLDDT 0.744 -> 0.823; pairwise epitope Jaccard 0.234 -> 0.250; hotspot coverage 0.301 -> 0.167.
- Native/Immunogen >=50% consensus Jaccard in the shared region: 0.000 -> 0.632 with MSA.
- MSA therefore materially changes the predicted interface. The 3-seed sensitivity run does not replace the 15-seed formal ensembles.

## Pose, glycosylation and geometry audit

- empty-MSA formal ensembles: after strict shared 27-160 antigen N/CA/C alignment, median antigen backbone RMSD=12.30 A and median antibody-pose C-alpha RMSD=42.24 A (all Native/Immunogen seed pairs).
- ColabFold-MSA formal ensembles: after strict shared 27-160 antigen N/CA/C alignment, median antigen backbone RMSD=3.50 A and median antibody-pose C-alpha RMSD=30.66 A (all Native/Immunogen seed pairs).
- ColabFold-MSA sensitivity ensembles: after strict shared 27-160 antigen N/CA/C alignment, median antigen backbone RMSD=1.97 A and median antibody-pose C-alpha RMSD=16.91 A (all Native/Immunogen seed pairs).
- The only sequence sequon is N74 (NVS). A matched single-NAG sensitivity ensemble lowered Native mean ipTM by 0.018 and complex pLDDT by 0.024; one Native model placed NAG only 1.51 A from the antibody, demonstrating steric incompatibility of that pose.
- Every unrelaxed Boltz model contains at least one severe <1.8 A nonbonded clash under the audit definition; use unrelaxed coordinates for hypothesis generation/ranking only.

## Relative affinity

- immunogen_blind: mean PRODIGY dG -14.69 kcal/mol (SD 1.38); relative ranking only.
- native_blind: mean PRODIGY dG -13.83 kcal/mol (SD 1.49); relative ranking only.
- immunogen_blind with formal MSA: mean PRODIGY dG -12.68 kcal/mol (SD 1.24); relative ranking only.
- native_blind with formal MSA: mean PRODIGY dG -12.35 kcal/mol (SD 0.97); relative ranking only.

## Restrained relaxation

- 5 representative complexes were minimized with Amber14/GBn2 while restraining backbone heavy atoms.
- immunogen_blind_blind_representative_seed1_relaxed: interchain clashes <1.8 A 7 -> 0; nonlocal intrachain clashes 1 -> 0; backbone coordinate RMSD 0.189 A.
- immunogen_blind_msa_consensus_representative_seed7_relaxed: interchain clashes <1.8 A 7 -> 0; nonlocal intrachain clashes 1 -> 0; backbone coordinate RMSD 0.108 A.
- native_blind_blind_representative_seed14_relaxed: interchain clashes <1.8 A 7 -> 0; nonlocal intrachain clashes 9 -> 0; backbone coordinate RMSD 0.135 A.
- native_blind_msa_consensus_representative_seed9_relaxed: interchain clashes <1.8 A 0 -> 0; nonlocal intrachain clashes 5 -> 0; backbone coordinate RMSD 0.191 A.
- native_blind_phase1_informed_representative_seed3_relaxed: interchain clashes <1.8 A 11 -> 0; nonlocal intrachain clashes 6 -> 0; backbone coordinate RMSD 0.191 A.

## Interpretation limits

- Formal results now include matched 15-seed empty-MSA and ColabFold-MSA ensembles; disagreement between conditions and seed-level multimodality remain important uncertainty.
- The N74-NAG test contains one core GlcNAc only, not a complete heterogeneous human N-glycan.
- Phase1-informed means post-ranking of the Native-Blind ensemble, not a directly restrained prediction.
- The available Phase 1B ranges came from the handoff record; the formal Phase 1B final candidate table was not present in phase2.
- Phase 1B is a prior, not ground truth; hotspot overlap is not accuracy or proof of an epitope.
- Protein-protein affinity values are relative modeling aids, not experimental KD.
- No reliable IL-24/receptor complex was present in phase2, so neutralization mechanism remains a hypothesis.
