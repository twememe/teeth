#!/usr/bin/env python3
"""Run the released entries and compare their real outputs with frozen references."""
import argparse,csv,datetime,hashlib,json,os,pathlib,subprocess,sys,time,traceback,zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with pathlib.Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 with pathlib.Path(p).open(newline='',encoding='utf-8-sig') as f:
  r=csv.DictReader(f);return r.fieldnames,list(r)
def require(ok,msg):
 if not ok:raise ValueError(msg)
def compare(a,b):
 require(a.is_file() and b.is_file(),'Missing comparison file: '+str(a)+' / '+str(b))
 if a.suffix=='.csv':require(rows(a)==rows(b),'CSV content mismatch: '+str(a))
 elif a.suffix=='.json':require(json.loads(a.read_text())==json.loads(b.read_text()),'JSON content mismatch: '+str(a))
 else:require(a.read_bytes()==b.read_bytes(),'File mismatch: '+str(a))
 return {'actual':str(a),'reference':str(b),'actual_sha256':sha(a),'reference_sha256':sha(b),'content_equal':True}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stage',choices=['cpu','base','esm','local'],required=True);ap.add_argument('--test-root',type=pathlib.Path,required=True);ap.add_argument('--zip',type=pathlib.Path,help='Optional submitted ZIP to bind to this extracted copy');a=ap.parse_args()
 require(a.test_root.is_absolute(),'--test-root must be absolute')
 T=a.test_root.resolve();T.mkdir(parents=True,exist_ok=True);L=T/'acceptance_logs';L.mkdir(exist_ok=True)
 recordfile=T/('acceptance_'+a.stage+'.json');require(not recordfile.exists(),'Stage already recorded; use a new test root for a rerun')
 env=dict(os.environ,PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg',OMP_NUM_THREADS='8',MKL_NUM_THREADS='8',HF_HOME=str(T/'hf_cache'),TORCH_HOME=str(T/'torch_cache'),XDG_CACHE_HOME=str(T/'xdg_cache'),TMPDIR=str(T))
 env.pop('PYTHONPATH',None);env['PATH']=str(pathlib.Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
 rec={'stage':a.stage,'cwd':str(ROOT),'python':sys.executable,'python_version':sys.version,'start':datetime.datetime.now(datetime.timezone.utc).isoformat(),'commands':[],'comparisons':[],'status':'RUNNING','policy':{'probability_atol':0,'probability_rtol':0,'basis':'Same RTX A6000, declared software, original eval/seed42/bf16 settings; fixed before execution','base_fingerprints':'exact','metadata_exclusions':[]},'environment_overrides':{k:env[k] for k in ['PYTHONNOUSERSITE','OMP_NUM_THREADS','MKL_NUM_THREADS','HF_HOME','TORCH_HOME','XDG_CACHE_HOME','TMPDIR']}}
 def save():recordfile.write_text(json.dumps(rec,ensure_ascii=False,indent=2)+'\n')
 def run(name,args,extra=None):
  c={'name':name,'command':[str(x) for x in args],'cwd':str(ROOT),'start':datetime.datetime.now(datetime.timezone.utc).isoformat(),'log':str(L/(name+'.log')),'environment_overrides':extra or {}};rec['commands'].append(c);save();start=time.monotonic();print(name,flush=True)
  with open(c['log'],'w') as f:p=subprocess.run(c['command'],cwd=ROOT,env=dict(env,**(extra or {})),stdout=f,stderr=subprocess.STDOUT)
  c.update(exit_code=p.returncode,elapsed_seconds=time.monotonic()-start,end=datetime.datetime.now(datetime.timezone.utc).isoformat());save();require(p.returncode==0,name+' failed; see '+c['log'])
 def cmp(p,q):rec['comparisons'].append(compare(p,q));save()
 def infer_check(out,reference):
  actual=out/'lora_validation_predictions.csv';cmp(actual,reference)
  _,pred=rows(actual);_,inp=rows(ROOT/'data/examples/esm2_sequences.csv');identity=json.loads((ROOT/'data/examples/esm2_manifest.json').read_text())
  require([r['row_id'] for r in pred]==identity['row_ids']==[r['row_id'] for r in inp],'Example row identity mismatch')
  required={r['row_id']:r for r in inp};found={}
  with (ROOT/'modules/inovate/data/validation.csv').open(newline='') as f:
   for r in csv.DictReader(f):
    if r['row_id'] in required:found[r['row_id']]=r
  require(set(found)==set(required),'Examples missing from saved validation split')
  for k,r in required.items():
   for col in ['model_sequence','label','source','record_id']:require(r[col]==found[k][col],'Sequence/identity mismatch: '+k+' '+col)
  _,refs=rows(reference);rec['inference']={'rows':len(pred),'row_ids':identity['row_ids'],'input_sha256':sha(ROOT/'data/examples/esm2_sequences.csv'),'max_probability_absolute_difference':max(abs(float(p['probability'])-float(q['probability'])) for p,q in zip(pred,refs)),'all_labels_equal':all(p['prediction']==q['prediction'] for p,q in zip(pred,refs)),'all_fields_equal':True,'raw_logits':'Not emitted by the formal evaluation CSV; raw_score is the source assay value, not a model logit.'}
 try:
  manifest={}
  for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
   digest,name=line.split(maxsplit=1);name=name.lstrip('*');require(sha(ROOT/name)==digest,'Release file SHA mismatch: '+name);manifest[name]=digest
  rec['release_files_verified']=len(manifest)
  if a.zip:
   zpath=a.zip.resolve();rec['zip']={'name':zpath.name,'path':str(zpath),'bytes':zpath.stat().st_size,'sha256':sha(zpath)}
   with zipfile.ZipFile(zpath) as z:
    candidates=[n for n in z.namelist() if n.endswith('/SHA256SUMS.txt')];require(len(candidates)==1,'ZIP manifest layout invalid');require(z.read(candidates[0])==(ROOT/'SHA256SUMS.txt').read_bytes(),'Extracted manifest differs from specified ZIP')
  rec['trained_weights']=[]
  for item in json.loads((ROOT/'logs/base_externalization/trained_weights_preserved.json').read_text())['files']:
   require(sha(ROOT/item['path'])==item['sha256'],'Project checkpoint changed: '+item['path']);rec['trained_weights'].append(item)
  save()
  if a.stage=='cpu':
   run('T1_export',['bash','run.sh','--module','export','--output',T/'export']);cmp(T/'export/results.csv',ROOT/'results/results.csv')
   from Bio.PDB import PDBParser
   from Bio.SeqUtils import seq1
   _,rr=rows(T/'export/results.csv');chains=[]
   for r in rr:
    p=(T/'export'/r['source_pdb']).resolve();require(p.is_relative_to(T/'export'),'Structure outside export directory');cmp(p,ROOT/'results'/r['source_pdb']);model=PDBParser(QUIET=True).get_structure(r['pose_id'],p)[0]
    for c,col in [('A','antigen_sequence'),('H','VH_sequence'),('L','VL_sequence')]:
     seq=''.join(seq1(x.resname) for x in model[c] if x.id[0]==' ');require(seq==r[col],'PDB chain differs: '+str(p)+' '+c)
    chains.append({'pose_id':r['pose_id'],'pdb':str(p),'sha256':sha(p),'A_H_L_sequences_match':True})
   rec['T1']={'candidate_count':len(rr),'unique_antibodies':len({r['antibody_id'] for r in rr}),'all_csv_fields_equal':True,'structures':chains};save()
   run('T2_phase1',[sys.executable,'run.py','--module','phase1','--output',T/'phase1'])
   p1=ROOT/'results/phase1'
   for q in sorted(p1.rglob('*.csv')):cmp(T/'phase1'/q.relative_to(p1),q)
   cmp(T/'phase1/integration/bepipred_discotope_correlations.json',p1/'integration/bepipred_discotope_correlations.json');cmp(T/'phase1/RUN_INFO.json',p1/'RUN_INFO.json')
   rec['T2']={'residues':len(rows(T/'phase1/integration/residue_feature_table_3d.csv')[1]),'families':len(rows(T/'phase1/candidate_regions.csv')[1]),'windows':{str(n):len(rows(T/f'phase1/window_{n}.csv')[1]) for n in [15,20,25]},'upstream':'Existing BepiPred/DiscoTope/SASA/DSSP/conservation evidence reused; no upstream prediction rerun'};save()
   run('T3_phase2',[sys.executable,'run.py','--module','phase2','--output',T/'phase2'])
   for name in ['contacts/complex_model_summary.csv','contacts/interface_contacts.csv','mechanism/6DF3_model_overlap.csv','mechanism/6DF3_receptor_contacts.json','mechanism/T198_project_mapping.json','mechanism/external_structure_analysis.md']:
    cmp(T/'phase2'/name,ROOT/'results/phase2'/name)
   for q in sorted((ROOT/'results/phase2/figures').glob('*.svg')):cmp(T/'phase2/figures'/q.name,q)
   rec['T3']={'models':len(rows(T/'phase2/contacts/complex_model_summary.csv')[1]),'residue_contact_pairs':len(rows(T/'phase2/contacts/interface_contacts.csv')[1]),'representatives':len(rows(T/'phase2/mechanism/6DF3_model_overlap.csv')[1]),'conditions':'Original PDB heavy-atom contact calculation, 4.5/5.0 Angstrom; 6DF3 protein-only IL-24-side interface and original numbering mapping; all result fields checked','no_Boltz_or_PRODIGY_rerun':True};save()
  elif a.stage=='base':
   cache=T/'official_cache';model=ROOT/'external_models/esm2_650m_hf';require(not cache.exists() or not any(cache.iterdir()),'T4 requires a new empty official cache');require(not model.exists(),'T4 requires no prepared model in extracted ZIP')
   rec['T4']={'empty_official_cache_before':True,'converted_model_absent_before':True};save()
   run('T4_prepare_base',[sys.executable,'tools/prepare_base.py','--cache-dir',cache])
   m=json.loads((ROOT/'models/base_manifest.json').read_text());rec['T4']['downloaded_sources']=[]
   for asset in m['source_files']:
    p=cache/asset['filename'];require(sha(p)==asset['sha256'] and p.stat().st_size==asset['bytes'],'Official source mismatch');rec['T4']['downloaded_sources'].append(asset)
   rec['T4']['preparation']=json.loads((model/'preparation_record.json').read_text());rec['T4']['reference_tensor_count']=len(json.loads((ROOT/'models/base_reference.json').read_text())['tensors'])
  else:
   run(a.stage+'_hardware',['nvidia-smi','--query-gpu=index,name,uuid,memory.used,utilization.gpu,driver_version','--format=csv'])
   gpu=subprocess.check_output(['nvidia-smi','--query-gpu=index,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True)
   available=[r.split(',')[0].strip() for r in gpu.splitlines() if int(r.split(',')[1])<512 and int(r.split(',')[2])==0];require(bool(available),'No idle GPU; retry later without interrupting other tasks');extra={'CUDA_VISIBLE_DEVICES':available[0]};rec['gpu_selected']=available[0]
   if a.stage=='esm':
    run('T5_esm2',[sys.executable,'run.py','--module','esm2','--output',T/'esm2'],extra);infer_check(T/'esm2',ROOT/'results/esm2/lora_validation_predictions.csv')
   else:
    run('T6_verify_local',[sys.executable,'tools/prepare_base.py','--converted-dir','./external_models/esm2_650m_hf','--verify-only'])
    extra.update(HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1');run('T6_esm2_local',[sys.executable,'run.py','--module','esm2','--model-dir','./external_models/esm2_650m_hf','--output',T/'esm2_local'],extra);infer_check(T/'esm2_local',ROOT/'results/esm2/lora_validation_predictions.csv');cmp(T/'esm2_local/lora_validation_predictions.csv',T/'esm2/lora_validation_predictions.csv')
  rec['status']='PASS'
 except Exception as e:
  rec['status']='FAIL';rec['error']=str(e);rec['traceback']=traceback.format_exc();raise
 finally:
  rec['end']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
 print(json.dumps({'status':rec['status'],'record':str(recordfile)},indent=2))
if __name__=='__main__':main()
