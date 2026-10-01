#!/usr/bin/env python3
"""Export the existing formal-MSA top-10 pose definition with unchanged sorting."""
import argparse,csv,json,shutil,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PACKAGE=ROOT.parents[1]

def fasta(p):
 d={};key=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):key=line[1:].split()[0];d[key]=''
  elif line.strip():d[key]+=line.strip()
 return d

def main():
 global ROOT
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--ranking',type=Path,default=ROOT/'results/14_quality/model_ranking_formal_msa.csv')
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--source-root',type=Path,default=ROOT)
 a=p.parse_args();ROOT=a.source_root.resolve();a.output.mkdir(parents=True,exist_ok=True)
 rows=list(csv.DictReader(a.ranking.open(encoding='utf-8')))
 required={'route','seed','confidence_score','iptm','complex_plddt','prodigy_dG_kcal_mol'}
 if not rows or not required.issubset(rows[0]):raise ValueError('Missing ranking rows/columns')
 for r in rows:
  for k in required-{'route','seed'}:
   if not math.isfinite(float(r[k])):raise ValueError(f'Invalid {k}')
 rows.sort(key=lambda r:(float(r['confidence_score']),float(r['iptm'])),reverse=True)
 seq=fasta(ROOT/'data/VH_VL_clean.fasta')
 vh=next(v for k,v in seq.items() if 'VH' in k);vl=next(v for k,v in seq.items() if 'VL' in k)
 version=json.loads((PACKAGE/'configs/release.json').read_text())['release_version']
 output=[];(a.output/'structures').mkdir(exist_ok=True)
 # Ten is the existing submission definition in the source predict.py, not a new contest rule.
 for i,r in enumerate(rows[:10],1):
  route,seed=r['route'],r['seed']
  source=ROOT/'results/12_complex_prediction/formal_msa'/f'{route}_seed{seed}'/f'boltz_results_{route}_msa/predictions/{route}_msa/{route}_msa_model_0.pdb'
  if not source.is_file():raise FileNotFoundError(source)
  relative=f'structures/{route}_seed{seed}.pdb';shutil.copy2(source,a.output/relative)
  native=route=='native_blind'
  antigen=next(iter(fasta(ROOT/'data'/('IL24_native_27_181.fasta' if native else 'IL24_immunogen_27_160.fasta')).values()))
  output.append({'candidate_id':f'IA6-13-8-{i:02d}','track':'赛道一，AI大分子与多肽药物设计','antibody_id':'IA6-13-8','pose_id':f'{route}_formal_msa_seed{seed}','route':route,'construct':'Native155 (27-181)' if native else 'Immunogen134 (27-160)','VH_sequence':vh,'VL_sequence':vl,'antigen_sequence':antigen,'seed':seed,'msa_condition':'formal_msa','sequence_or_structure':relative,'source_pdb':relative,'structure_status':'unrelaxed','coordinate_unit':'angstrom','model_version':'Boltz-1 weights; boltz 2.2.1','run_version':version,'confidence_score':r['confidence_score'],'ipTM':r['iptm'],'complex_pLDDT':r['complex_plddt'],'PRODIGY_dG_kcal_mol':r['prodigy_dG_kcal_mol'],'hotspot_coverage':r['hotspot_coverage'],'ranking_basis':f"confidence_rank={r['confidence_rank']}; consensus_rank={r['consensus_rank']}",'notes':'Same antibody molecular identity across rows; distinct computational poses. Confidence/ipTM/pLDDT/coverage are dimensionless (Boltz complex_pLDDT is 0-1). PRODIGY is predicted dG, not measured affinity. Raw coordinates may contain clashes; see Model Card.'})
 with (a.output/'results.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
 print('Exported',len(output),'existing poses to',a.output/'results.csv')

if __name__=='__main__':main()
