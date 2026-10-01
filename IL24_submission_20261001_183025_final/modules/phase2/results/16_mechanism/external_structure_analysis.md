# Phase 2 external-structure analysis

Corrected computation (logic closeout 2026-09-30).  The previous version of this file and of the
derived tables used the wrong object: it admitted HETATM records into the contact map, applied the
construct-local +26 convention to 6DF3 chain C, and intersected IL-24 epitope numbers with
receptor-chain numbers.  Those outputs remain in the original engineering archive and are excluded from this release.

## What was computed

- Protein contacts only: amino-acid `ATOM` residues, heavy atoms, 5.0 A cutoff.  HOH, NAG
  and GOL are excluded.
- Chain identity read from the structure: **chain C = human IL-24** (UNP Q13007 52-206),
  **chain L = IL-22RA1** (Q8N6P7 24-228), **chain H = IL-20RB** (Q6UXL0 35-224).  In this external
  complex L and H are *receptors*.
- 6DF3 chain C substitutions relative to canonical human IL-24: N85Q, N99Q, Y124H, N126Q
  (the construct is annotated MUTATION: YES).
- **R** = the IL-24 side of the interface: 48 residues in 6DF3 author numbering, all of
  which map to project numbering (48 residues).  Receptor-side contacts
  (49 residues) are reported
  separately and are never used for the epitope overlap.
- Cross-species transfer: 6DF3 chain C aligned to the project reference 27-181.  The alignment is
  gap-free over all 155 columns with 106/155 identical positions
  (0.684), so project = 6DF3_auth - 25.  Of the 48 mapped contact
  residues, 34 have the same side chain in human and mouse.
- **E** = the model's IL-24 epitope, mapped to project numbering by locating its antigen-chain
  sequence inside the reference (no hard-coded offset).

## Result

`6DF3_model_overlap.csv` reports Jaccard = |E n R| / |E u R| for each shortlisted representative.
This is a genuine same-reference residue-set comparison.  It remains a structural hypothesis about
where the modelled antibody binds relative to the receptor footprint of a *different species* and
an *engineered* IL-24 construct; it is not proof of receptor competition, and the cross-species and
construct-mutation caveats above apply to every row.

- T198 maps to project residue 173 (verified by the
  alignment above, not by an assumed offset).
- ELISA data were preserved unchanged and summarized with antibody-minus-isotype inhibition in
  `ELISA_external_consistency.json`.
- 6DF3 does not establish the IL-20R1/IL-20R2 complex.
