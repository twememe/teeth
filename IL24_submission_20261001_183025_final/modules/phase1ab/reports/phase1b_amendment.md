# Phase 1B Protocol Amendment

## Status

This is a **post-freeze, user-authorized restoration of a preregistered missing evidence branch before historical unblinding**.

The original DiscoTope-3.0 branch was preregistered but remained missing in Phase 1A after external dependency/DNS retrieval failure exhausted the original installation ceiling. The user has now explicitly authorized a separate Phase 1B attempt to restore only that branch, then rerun integration and candidate extraction in new output directories.

## Boundaries

- `historical_region_still_blinded: true`
- `phase1a_outputs_mutable: false`
- `no_new_predictor_substitution: true`
- Step 0/1/2/3/4/6 are reused read-only and are not rerun.
- Phase 1A scoring, candidates, figures, and report are not overwritten.
- Step 9 is not authorized and will not be performed.

## New analysis

1. Create isolated `il24-dt3-phase1b` environment at official DiscoTope-3.0 commit `35d9f2e55f97eaba2a7acefbc394db58fb9670bc` where installable.
2. Run official AlphaFold structure mode and require exact 181/181 residue mapping.
3. Compute Phase 1B grouped evidence score using `E=(B+D)/2` and `Final=(E+R+C)/3`.
4. Rerun 15/20/25-aa candidate extraction, evidence-group ablation, model sensitivity, and Phase 1A-vs-Phase 1B comparison without historical data.
