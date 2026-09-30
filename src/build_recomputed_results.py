#!/usr/bin/env python3
"""Build an isolated standardized result table from a completed full run."""
import csv, os, shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ["PHASE2_OUTPUT_ROOT"]).resolve()
ranking=OUT/"14_quality/model_ranking_formal_msa.csv"
rows=list(csv.DictReader(ranking.open(encoding="utf-8")))
rows.sort(key=lambda r:(float(r["confidence_score"]),float(r["iptm"])),reverse=True)
candidates=OUT/"candidates"; candidates.mkdir(parents=True,exist_ok=True)
result=[]
for i,r in enumerate(rows[:10],1):
    route,seed=r["route"],r["seed"]
    source=next((OUT/"12_complex_prediction/formal_msa"/f"{route}_seed{seed}").glob("**/*_model_0.pdb"))
    dest=candidates/f"{route}_seed{seed}.pdb"; shutil.copy2(source,dest)
    result.append({"candidate_id":f"IA6-13-8-{i:02d}","track":"赛道一：AI大分子与多肽药物设计","task":"IL-24 antibody structure prediction","route":route,"seed":seed,"msa_condition":"formal_msa","sequence_or_structure":str(dest.relative_to(OUT)),"source_pdb":str(dest.relative_to(OUT)),"structure_status":"unrelaxed","relaxation_status":"not_relaxed","clash_qc":"raw model may contain severe clashes","model_version":f"Boltz-1 2.2.1; formal MSA; seed {seed}","ipTM":r["iptm"],"complex_pLDDT":r["complex_plddt"],"PRODIGY_dG_kcal_mol":r["prodigy_dG_kcal_mol"],"hotspot_coverage":r["hotspot_coverage"],"ranking_basis":f"confidence_rank={r['confidence_rank']}; consensus_rank={r['consensus_rank']}","notes":"Computational ranking only; PRODIGY is relative modeling value, not experimental KD."})
with (OUT/"results.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(result[0])); w.writeheader(); w.writerows(result)
print(f"Wrote {OUT/'results.csv'} with {len(result)} candidates")
