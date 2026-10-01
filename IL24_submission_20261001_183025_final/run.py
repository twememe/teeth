#!/usr/bin/env python3
"""Release entry: real export, evidence fusion, structural contacts and saved-model inference."""
import argparse,csv,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def call(args,cwd=None):
 print('+', ' '.join(map(str,args)),flush=True)
 subprocess.run([str(x) for x in args],cwd=cwd,check=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg'))

def phase1(inp,out):
 sys.path.insert(0,str(ROOT/'modules/phase1ab/scripts'))
 import pandas as pd
 import integration_phase1b as integration
 import candidates_phase1b as candidates
 from types import SimpleNamespace
 integration.run(SimpleNamespace(reference=inp/'inputs/reference_181.fasta',bepipred=inp/'results/03_bepipred/bepipred_residue_scores.csv',discotope=inp/'results/05b_discotope_restored/discotope_residue_scores.csv',plddt=inp/'results/02_structure_qc/residue_plddt.csv',sasa=inp/'results/04_surface/residue_sasa.csv',conservation=inp/'results/06_conservation/residue_conservation.csv',secondary=inp/'results/04_surface/residue_secondary_structure.csv',phase1a_table=inp/'results/07_integration/residue_feature_table.csv',output_dir=out/'integration',report=out/'integration.md'))
 table=pd.read_csv(out/'integration/residue_feature_table_3d.csv')
 settings=candidates.score_settings_phase1b(table)
 windows,families,scales=candidates._family_collections(settings)
 top=families['Full'][:5]
 stability=candidates._stability_matrix(top,families)
 scale_matches=candidates._scale_matches(top,scales)
 four=candidates._matching_table(top,families['4way-sensitivity'],setting='4way-sensitivity')
 ranking=candidates._candidate_table(table,top,stability,scale_matches,four)
 ranking.to_csv(out/'candidate_regions.csv',index=False)
 stability.to_csv(out/'candidate_stability.csv',index=False)
 for scale,rows in scales.items():candidates._window_frame(table,rows).to_csv(out/f'window_{scale}.csv',index=False)
 (out/'RUN_INFO.json').write_text(json.dumps({'mode':'real_fusion_window_analysis_and_ranking','upstream_reused':['official BepiPred residue outputs','official DiscoTope residue outputs','FreeSASA/DSSP','MAFFT conservation'],'formula':'(((B+D)/2)+R+C)/3','bootstrap':'Historical uncertainty is supplied in the module reference results; no new bootstrap was run.'},indent=2)+'\n')
 print('Phase 1:',len(table),'residues;',len(ranking),'ranked region families')

def phase2(inp,out):
 module=ROOT/'modules/phase2'
 call([sys.executable,module/'scripts/analyze_boltz_complexes.py','--set','formal_msa','--input',inp/'results/12_complex_prediction/formal_msa','--output',out/'contacts'])
 call([sys.executable,module/'scripts/run_phase2_external_analysis.py','--reference',inp/'data/IL24_full_1_181.fasta','--complex',inp/'data/external_materials/6DF3.pdb','--structures',inp/'results/14_quality/representatives_msa','--orthologs',ROOT/'modules/phase1ab/results/06_conservation/orthologs.fasta','--out',out/'mechanism'])
 call([sys.executable,module/'scripts/make_phase2_visualizations.py','--project-root',inp,'--mechanism',out/'mechanism','--out',out/'figures'])

def esm2(inp,out,config,base,checkpoint):
 if not (base/'config.json').is_file():raise FileNotFoundError('Prepare the frozen base first: python tools/prepare_base.py; or pass --model-dir with a verified converted directory')
 sys.path.insert(0,str(ROOT/'modules/inovate/scripts'))
 import torch,yaml
 import train
 cfg=yaml.safe_load(config.read_text())
 cfg.update(model_dir=str(base.resolve()),output_dir=str(out.resolve()),validation_csv=str(inp.resolve()),gradient_checkpointing=False)
 (out/'run_config.json').write_text(json.dumps(cfg,indent=2)+'\n')
 call([sys.executable,ROOT/'modules/inovate/scripts/train.py','--mode','evaluate','--config',out/'run_config.json','--checkpoint',checkpoint,'--split','validation','--max-seconds','600'])

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--module',choices=['export','phase1','phase2','esm2'],default='export')
 p.add_argument('--input',type=Path,help='export: ranking CSV; phase1/phase2: module-shaped input directory; esm2: CSV containing model_sequence,label,source')
 p.add_argument('--source-root',type=Path,help='export: matching Phase2 project root for a custom ranking CSV and its actual structures/data')
 p.add_argument('--output',type=Path,default=Path('run_output'),help='New, empty output directory; reference results are never overwritten')
 p.add_argument('--config',type=Path,default=ROOT/'configs/esm2_inference.json')
 p.add_argument('--model-dir',type=Path,default=ROOT/'external_models/esm2_650m_hf')
 p.add_argument('--checkpoint',type=Path,default=ROOT/'modules/inovate/output/esm2_650m_lora_best.pt')
 a=p.parse_args();out=a.output.resolve()
 if out.exists() and any(out.iterdir()):p.error('Output directory must be empty: '+str(out))
 for protected in [ROOT/'modules',ROOT/'models',ROOT/'results',ROOT/'data']:
  if out==protected or protected in out.parents:p.error('Choose an output directory outside reference material')
 out.mkdir(parents=True,exist_ok=True)
 if a.module=='export':
  source=(a.source_root or ROOT/'modules/phase2').resolve()
  rank=(a.input or source/'results/14_quality/model_ranking_formal_msa.csv').resolve()
  if source not in rank.parents:p.error('Custom ranking requires --source-root for its matching data and structures')
  call([sys.executable,ROOT/'modules/phase2/predict.py','--ranking',rank,'--source-root',source,'--output',out])
 elif a.module=='phase1':phase1((a.input or ROOT/'modules/phase1ab').resolve(),out)
 elif a.module=='phase2':phase2((a.input or ROOT/'modules/phase2').resolve(),out)
 else:esm2((a.input or ROOT/'data/examples/esm2_sequences.csv').resolve(),out,a.config.resolve(),a.model_dir,a.checkpoint)

if __name__=='__main__':main()
