#!/usr/bin/env python3
"""Create route-level ensemble summaries and a concise Phase 2 report."""
import csv, statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"results"/"13_interface_analysis"/"formal"/"complex_model_summary.csv"
OUT=ROOT/"results"/"17_delivery"; OUT.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader(SRC.open()))
aff_rows=list(csv.DictReader(open(ROOT/"results"/"15_affinity"/"prodigy_formal.csv")))
msa_rows=list(csv.DictReader(open(ROOT/"results"/"14_quality"/"msa_sensitivity_summary.csv")))
formal_msa_rows=list(csv.DictReader(open(ROOT/"results"/"14_quality"/"formal_msa_comparison.csv")))
glyco_rows=list(csv.DictReader(open(ROOT/"results"/"14_quality"/"glyco_sensitivity_summary.csv")))
pose_rows=list(csv.DictReader(open(ROOT/"results"/"14_quality"/"native_immunogen_pose_summary.csv")))
aff_msa_rows=list(csv.DictReader(open(ROOT/"results"/"15_affinity"/"prodigy_formal_msa.csv")))
by=defaultdict(list)
for r in rows: by[r["route"]].append(r)
metrics=["confidence_score","ptm","iptm","complex_plddt","epitope_residues_4p5","paratope_residues_4p5","cdr_paratope_residues_4p5","hotspot_coverage","hotspot_fraction_in_interface"]
agg=[]
for route,rs in by.items():
    x={"route":route,"n_models":len(rs)}
    for m in metrics:
        vals=[float(r[m]) for r in rs]; x[m+"_mean"]=statistics.mean(vals); x[m+"_sd"]=statistics.stdev(vals) if len(vals)>1 else 0
    counts=Counter(n for r in rs for n in r["epitope_project_residues"].split(";") if n)
    x["epitope_consensus_ge50pct"]=";".join(sorted((n for n,c in counts.items() if c>=len(rs)/2),key=int))
    x["epitope_consensus_ge30pct"]=";".join(sorted((n for n,c in counts.items() if c>=len(rs)*.3),key=int))
    agg.append(x)
with open(OUT/"route_ensemble_summary.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=agg[0].keys());w.writeheader();w.writerows(agg)

report=["# Phase 2 computational summary","",f"Formal models analyzed: {len(rows)} empty-MSA + {sum(int(r['n_models']) for r in formal_msa_rows if r['condition']=='formal_msa')} ColabFold-MSA; plus 6 N74-NAG sensitivity models.","","## Empty-MSA formal route summaries",""]
for x in agg:
    report += [f"### {x['route']}","",f"- n={x['n_models']}; mean ipTM={x['iptm_mean']:.3f} (SD {x['iptm_sd']:.3f})",f"- mean complex pLDDT={x['complex_plddt_mean']:.3f}",f"- mean Phase-1 hotspot coverage={x['hotspot_coverage_mean']:.3f}",f"- >=50% epitope consensus residues (project numbering): {x['epitope_consensus_ge50pct'] or 'none'}",""]
if len(agg)==2:
    sets=[set(x["epitope_consensus_ge50pct"].split(";"))-{""} for x in agg]
    jac=len(sets[0]&sets[1])/len(sets[0]|sets[1]) if sets[0]|sets[1] else 0
    report += ["## Native versus immunogen", "", f"- Jaccard index of >=50% epitope consensus sets: {jac:.3f}", "- Comparison uses project numbering and the shared 27-160 region; residue 163 in Native is outside the immunogen construct.", ""]
report += ["## Formal 15-seed MSA ensembles", ""]
for route in ("native_blind","immunogen_blind"):
    empty=next(r for r in formal_msa_rows if r["condition"]=="formal" and r["route"]==route)
    msa=next(r for r in formal_msa_rows if r["condition"]=="formal_msa" and r["route"]==route)
    report += [f"- {route}: mean ipTM {float(empty['iptm_mean']):.3f} -> {float(msa['iptm_mean']):.3f}; complex pLDDT {float(empty['complex_plddt_mean']):.3f} -> {float(msa['complex_plddt_mean']):.3f}; mean pairwise epitope Jaccard {float(empty['pairwise_epitope_jaccard_mean']):.3f} -> {float(msa['pairwise_epitope_jaccard_mean']):.3f}; hotspot coverage {float(empty['hotspot_coverage_mean']):.3f} -> {float(msa['hotspot_coverage_mean']):.3f}."]
msa_cross=next(r for r in formal_msa_rows if r["condition"]=="formal_msa")["native_immunogen_consensus_jaccard"]
report += [f"- Native/Immunogen >=50% consensus Jaccard in the shared region is {float(msa_cross):.3f} with formal MSA, versus {jac:.3f} with empty MSA.", "- The formal-MSA consensus is concentrated at project residues 132-150, but seed-level pose heterogeneity remains substantial.", ""]
report += ["## Matched MSA sensitivity analysis", "", "Three matched seeds per route used 50 sampling steps for both empty-MSA and ColabFold-MSA conditions.", ""]
for route in ("native_blind","immunogen_blind"):
    empty=next(r for r in msa_rows if r["route"]==route and r["condition"]=="empty_msa")
    msa=next(r for r in msa_rows if r["route"]==route and r["condition"]=="colabfold_msa")
    report += [f"- {route}: ipTM {float(empty['iptm_mean']):.3f} -> {float(msa['iptm_mean']):.3f}; complex pLDDT {float(empty['complex_plddt_mean']):.3f} -> {float(msa['complex_plddt_mean']):.3f}; pairwise epitope Jaccard {float(empty['pairwise_epitope_jaccard']):.3f} -> {float(msa['pairwise_epitope_jaccard']):.3f}; hotspot coverage {float(empty['hotspot_coverage_mean']):.3f} -> {float(msa['hotspot_coverage_mean']):.3f}."]
empty_cross=next(r for r in msa_rows if r["condition"]=="empty_msa")["native_immunogen_consensus_jaccard"]
msa_cross=next(r for r in msa_rows if r["condition"]=="colabfold_msa")["native_immunogen_consensus_jaccard"]
report += [f"- Native/Immunogen >=50% consensus Jaccard in the shared region: {float(empty_cross):.3f} -> {float(msa_cross):.3f} with MSA.", "- MSA therefore materially changes the predicted interface. The 3-seed sensitivity run does not replace the 15-seed formal ensembles.", ""]
report += ["## Pose, glycosylation and geometry audit", ""]
for condition,label in (("formal","empty-MSA formal ensembles"),("formal_msa","ColabFold-MSA formal ensembles"),("msa_sensitivity","ColabFold-MSA sensitivity ensembles")):
    antigen=next(r for r in pose_rows if r["condition"]==condition and r["metric"]=="shared_antigen_backbone_rmsd_A")
    antibody=next(r for r in pose_rows if r["condition"]==condition and r["metric"]=="antibody_pose_ca_rmsd_A")
    report += [f"- {label}: after strict shared 27-160 antigen N/CA/C alignment, median antigen backbone RMSD={float(antigen['median']):.2f} A and median antibody-pose C-alpha RMSD={float(antibody['median']):.2f} A (all Native/Immunogen seed pairs)."]
report += ["- The only sequence sequon is N74 (NVS). A matched single-NAG sensitivity ensemble lowered Native mean ipTM by 0.018 and complex pLDDT by 0.024; one Native model placed NAG only 1.51 A from the antibody, demonstrating steric incompatibility of that pose.", "- Every unrelaxed Boltz model contains at least one severe <1.8 A nonbonded clash under the audit definition; use unrelaxed coordinates for hypothesis generation/ranking only.", ""]
report += ["## Relative affinity",""]
for route in sorted({r["route"] for r in aff_rows}):
    vals=[float(r["prodigy_dG_kcal_mol"]) for r in aff_rows if r["route"]==route]
    report += [f"- {route}: mean PRODIGY dG {statistics.mean(vals):.2f} kcal/mol (SD {statistics.stdev(vals):.2f}); relative ranking only."]
for route in sorted({r["route"] for r in aff_msa_rows}):
    vals=[float(r["prodigy_dG_kcal_mol"]) for r in aff_msa_rows if r["route"]==route]
    report += [f"- {route} with formal MSA: mean PRODIGY dG {statistics.mean(vals):.2f} kcal/mol (SD {statistics.stdev(vals):.2f}); relative ranking only."]
report += [""]
relaxed_path=ROOT/"results"/"14_quality"/"relaxed_geometry_qc.csv"
if relaxed_path.exists():
    relaxed=list(csv.DictReader(relaxed_path.open()))
    report += ["## Restrained relaxation", "", f"- {len(relaxed)} representative complexes were minimized with Amber14/GBn2 while restraining backbone heavy atoms."]
    for r in relaxed:
        report += [f"- {r['model']}: interchain clashes <1.8 A {r['before_interchain_clashes_lt1p8']} -> {r['after_interchain_clashes_lt1p8']}; nonlocal intrachain clashes {r['before_nonlocal_intrachain_clashes_lt1p8']} -> {r['after_nonlocal_intrachain_clashes_lt1p8']}; backbone coordinate RMSD {float(r['backbone_coordinate_rmsd_A']):.3f} A."]
    report += [""]
report += ["## Interpretation limits","","- Formal results now include matched 15-seed empty-MSA and ColabFold-MSA ensembles; disagreement between conditions and seed-level multimodality remain important uncertainty.","- The N74-NAG test contains one core GlcNAc only, not a complete heterogeneous human N-glycan.","- Phase1-informed means post-ranking of the Native-Blind ensemble, not a directly restrained prediction.","- The available Phase 1B ranges came from the handoff record; the formal Phase 1B final candidate table was not present in phase2.","- Phase 1B is a prior, not ground truth; hotspot overlap is not accuracy or proof of an epitope.","- Protein-protein affinity values are relative modeling aids, not experimental KD.","- No reliable IL-24/receptor complex was present in phase2, so neutralization mechanism remains a hypothesis.",""]
(OUT/"Phase2_summary.md").write_text("\n".join(report))
print(OUT)
