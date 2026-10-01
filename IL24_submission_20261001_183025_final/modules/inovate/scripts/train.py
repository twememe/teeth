"""Deadline-aware ESM-2 LoRA experiment. Test data is read only with --split test."""
import argparse, contextlib, csv, hashlib, json, math, os, random, sys, time, traceback
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.distributed as dist
from scipy.stats import spearmanr
from sklearn.metrics import f1_score, matthews_corrcoef
import yaml

ROOT=Path(__file__).resolve().parents[1]
RANK=0; WORLD=1; DEVICE=torch.device('cpu'); OUT=ROOT/'output'


def setup():
 global RANK,WORLD,DEVICE
 for k,v in json.loads((ROOT/'environment.json').read_text()).items():os.environ.setdefault(k,v)
 WORLD=int(os.getenv('WORLD_SIZE','1'));RANK=int(os.getenv('RANK','0'))
 if torch.cuda.is_available():
  local=int(os.getenv('LOCAL_RANK','0'));torch.cuda.set_device(local);DEVICE=torch.device('cuda',local)
 if WORLD>1:dist.init_process_group('nccl' if DEVICE.type=='cuda' else 'gloo',timeout=timedelta(minutes=15))


def log(message):
 text=f'{datetime.now(timezone.utc).isoformat()} rank={RANK} {message}'
 print(text,flush=True)
 OUT.mkdir(parents=True,exist_ok=True)
 with (OUT/f'run_rank{RANK}.log').open('a') as f:f.write(text+'\n')
 if RANK==0:
  with (OUT/'progress.txt').open('a') as f:f.write(text+'\n')


def seed_all(seed):
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 if torch.cuda.is_available():torch.cuda.manual_seed_all(seed)


def barrier():
 if WORLD>1:dist.barrier()


def deadline_reached(deadline,collective=False):
 stop=time.time()>=deadline
 if collective and WORLD>1:
  t=torch.tensor(int(stop),device=DEVICE);dist.all_reduce(t,op=dist.ReduceOp.MAX);stop=bool(t.item())
 return stop


def gather(objects):
 if WORLD==1:return objects
 chunks=[None]*WORLD;dist.all_gather_object(chunks,objects)
 return [v for chunk in chunks for v in chunk]


def metrics(y,p):
 y=np.asarray(y,dtype=int);p=np.asarray(p,dtype=float)
 if len(y)==0:return {'n':0,'mcc':None,'f1':None,'spearman_probability_label':None}
 rho=float(spearmanr(p,y).statistic) if len(set(y))>1 and len(set(p))>1 else None
 return {'n':len(y),'positive_fraction':float(y.mean()),'mcc':float(matthews_corrcoef(y,p>=.5)),
         'f1':float(f1_score(y,p>=.5,zero_division=0)),'spearman_probability_label':rho,'threshold':.5}


def grouped_metrics(frame,probabilities):
 result={'overall':metrics(frame.label,probabilities),'by_source':{},'by_source_partner':{},'within_source_raw_score_spearman':{}}
 p=np.asarray(probabilities)
 for source,idx in frame.groupby('source',sort=True).indices.items():
  result['by_source'][str(source)]=metrics(frame.label.iloc[idx],p[idx])
  if 'raw_score' in frame:
   raw=pd.to_numeric(frame.raw_score.iloc[idx],errors='coerce').to_numpy();valid=np.isfinite(raw)
   rho=float(spearmanr(raw[valid],p[idx][valid]).statistic) if valid.sum()>2 and len(np.unique(raw[valid]))>1 and len(np.unique(p[idx][valid]))>1 else None
   result['within_source_raw_score_spearman'][str(source)]={'n':int(valid.sum()),'spearman':rho,'interpretation':'Numerical raw-score direction; heterogeneous assays may differ.'}
 if 'mutated_partner' in frame:
  for (source,partner),idx in frame.groupby(['source','mutated_partner'],dropna=False).indices.items():
   result['by_source_partner'][str(source)+'/'+str(partner)]=metrics(frame.label.iloc[idx],p[idx])
 return result


def load_frame(path):
 frame=pd.read_csv(path,keep_default_na=False,low_memory=False)
 required={'model_sequence','label','source'}
 if not required.issubset(frame):raise ValueError(f'Missing required columns: {required-set(frame.columns)}')
 sequences=frame.model_sequence.astype(str)
 invalid=(sequences.str.len()<1)|(sequences.str.len()>1022)|~sequences.str.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZOU]+')
 if invalid.any():raise ValueError(f'{path}: {int(invalid.sum())} invalid or >1022-residue sequences; never truncate silently')
 frame['label']=pd.to_numeric(frame.label)
 if not set(frame.label).issubset({0,1}):raise ValueError('Labels must be binary')
 frame['sample_index']=np.arange(len(frame))
 if 'row_id' not in frame:frame['row_id']=[f'{Path(path).name}:{i}' for i in range(len(frame))]
 return frame


def split_indices(n,rank,world):return list(range(rank,n,world))


def bucket_batches(lengths,batch_size,rank,world,seed,training):
 rng=np.random.default_rng(seed);indices=np.arange(len(lengths));lengths=np.asarray(lengths)
 if training:
  rng.shuffle(indices);window=batch_size*world*50
  indices=np.concatenate([a[np.argsort(lengths[a],kind='stable')] for a in [indices[i:i+window] for i in range(0,len(indices),window)]]) if len(indices) else indices
  width=batch_size*world;chunks=[indices[i:i+width] for i in range(0,len(indices)-width+1,width)]
  rng.shuffle(chunks)
  return [x[rank*batch_size:(rank+1)*batch_size].tolist() for x in chunks]
 indices=np.asarray(split_indices(len(indices),rank,world));indices=indices[np.argsort(lengths[indices],kind='stable')]
 return [indices[i:i+batch_size].tolist() for i in range(0,len(indices),batch_size)]


class Tokens:
 def __init__(self,frame,model_dir):
  from transformers import AutoTokenizer
  tokenizer=AutoTokenizer.from_pretrained(model_dir,local_files_only=True)
  self.pad=tokenizer.pad_token_id
  self.tokens=tokenizer(frame.model_sequence.tolist(),add_special_tokens=True,truncation=False,padding=False)['input_ids']
  if any(len(ids)!=len(seq)+2 for ids,seq in zip(self.tokens,frame.model_sequence)):raise ValueError('Tokenizer did not produce one token per residue')
  self.labels=frame.label.to_numpy(dtype=np.float32);self.lengths=list(map(len,self.tokens))
 def batch(self,indices):
  x=torch.full((len(indices),max(len(self.tokens[i]) for i in indices)),self.pad,dtype=torch.long)
  for j,i in enumerate(indices):x[j,:len(self.tokens[i])]=torch.tensor(self.tokens[i])
  x=x.to(DEVICE);return x,x.ne(self.pad).long(),torch.tensor(self.labels[indices],device=DEVICE)


def amp():return torch.autocast('cuda',dtype=torch.bfloat16) if DEVICE.type=='cuda' else contextlib.nullcontext()


def atomic_json(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False));tmp.replace(path)


def atomic_save(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');torch.save(obj,tmp);tmp.replace(path)


def sha_file(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()


def source_identity(config):
 # Validate the loaded runtime tensors/config/tokenizer, independent of serialization bytes.
 sys.path.insert(0,str(ROOT.parents[1]/'tools'))
 from prepare_base import verify_converted
 return verify_converted(config['model_dir'])


def cache_matches(metadata,data,model):
 return metadata.get('data')==data and metadata.get('model')==model


def require_model_identity(actual,expected,context):
 sys.path.insert(0,str(ROOT.parents[1]/'tools'))
 from prepare_base import normalize_identity
 if normalize_identity(actual)!=normalize_identity(expected):
  raise ValueError(f'{context}: model differs from the frozen reference')


def data_identity(path):return {'path':str(path),'sha256':sha_file(path)}


def write_predictions(name,frame,p,extra=None):
 if RANK:return
 cols=[c for c in ['sample_index','row_id','record_id','source','assay','mutated_partner','label','raw_score','label_semantics'] if c in frame]
 pred=frame[cols].copy();pred['probability']=p;pred['prediction']=(np.asarray(p)>=.5).astype(int)
 pred.to_csv(OUT/f'{name}_validation_predictions.csv',index=False)
 report=grouped_metrics(frame,p);report.update(extra or {})
 atomic_json(OUT/f'{name}_metrics.json',report)
 report_prefix='smoke_' if name.startswith('smoke_') else ''
 summary_path=OUT/(report_prefix+'metrics.json');summary=json.loads(summary_path.read_text()) if summary_path.exists() else {}
 summary[name]=report
 summary['theoretical_balanced_random']={'f1':.5,'mcc':0.,'note':'Population expectation for independent 50/50 labels and predictions; not measured results.'}
 summary['label_interpretation']='Numerical high score within each source/assay; biological directions are not assumed equivalent across sources.'
 atomic_json(summary_path,summary)
 rows=[]
 for key,value in summary.items():
  if isinstance(value,dict) and 'overall' in value:rows.append({'model':key,**value['overall']})
 pd.DataFrame(rows).to_csv(OUT/(report_prefix+'comparison.csv'),index=False)
 lines=['# Measured validation comparison','','All rows use the same validation partition unless explicitly marked test. Labels represent numerical high scores, not a unified biological positive direction.','','| Model | N | MCC | F1 | Spearman (probability, label) |','|---|---:|---:|---:|---:|']
 for row in rows:lines.append('| '+ ' | '.join(str(row.get(k)) for k in ['model','n','mcc','f1','spearman_probability_label'])+' |')
 lines+=['','Theoretical balanced random reference: F1 = 0.5, MCC = 0 (population expectations, not measurements).','See each metrics JSON for source and mutated-partner subgroups.']
 (OUT/(report_prefix+'comparison.md')).write_text('\n'.join(lines)+'\n')
 if name in {'lora','smoke_lora'}:pred.to_csv(OUT/('smoke_validation_predictions.csv' if name.startswith('smoke_') else 'validation_predictions.csv'),index=False)


def random_baseline(frame,config,prefix=''):
 if RANK==0:write_predictions(prefix+'random',frame,np.random.default_rng(config['seed']).random(len(frame)),{'seed':config['seed'],'method':'Independent uniform [0,1) probabilities, threshold 0.5; measured on identical validation rows.'})


def stitch_parts(paths,n):
 output=None;seen=np.zeros(n,dtype=bool)
 for path in paths:
  with np.load(path) as part:
   idx=part['indices'];features=part['features']
   if output is None:output=np.empty((n,features.shape[1]),dtype=np.float32)
   if len(idx)!=len(features) or len(np.unique(idx))!=len(idx) or np.any(idx<0) or np.any(idx>=n) or seen[idx].any():raise ValueError('Duplicate or invalid extraction indices')
   output[idx]=features;seen[idx]=True
 if not seen.all():raise ValueError(f'Incomplete extraction: {seen.sum()}/{n}')
 return output


def extraction(config,deadline,smoke=False):
 from modeling import make_model
 prefix='smoke_' if smoke else '';model_identity=source_identity(config);model=make_model(config,lora=False).to(DEVICE);model.eval()
 for split in ['train','validation']:
  path=config[(prefix+split+'_csv')];frame=load_frame(path);tokens=Tokens(frame,config['model_dir']);identity=data_identity(path)
  cache=OUT/f'{prefix}{split}_embeddings.npy';meta=OUT/f'{prefix}{split}_embeddings.json'
  if cache.exists() and meta.exists() and cache_matches(json.loads(meta.read_text()),identity,model_identity):
   log(f'Using verified embedding cache {cache}');continue
  collected=[];indices=[]
  batches=bucket_batches(tokens.lengths,config['extraction_batch_size'],RANK,WORLD,config['seed'],False)
  with torch.inference_mode():
   for step,idx in enumerate(batches):
    if deadline_reached(deadline):break
    x,mask,_=tokens.batch(idx)
    with amp():features=model(x,mask,features_only=True)
    collected.append(features.float().cpu().numpy());indices.extend(idx)
    if step%20==0:log(f'extract {split} batch={step+1}/{len(batches)} local_rows={len(indices)}')
  arr=np.concatenate(collected) if collected else np.empty((0,model.head.in_features),dtype=np.float32)
  part=OUT/f'{prefix}{split}_embedding_part_rank{RANK}.npz';np.savez(part,indices=np.asarray(indices,dtype=int),features=arr)
  complete=len(indices)==len(split_indices(len(frame),RANK,WORLD));statuses=gather([complete])
  if not all(statuses):log(f'Extraction stopped at deadline; {split} cache incomplete and will not be consumed');return False
  barrier()
  if RANK==0:
   array=stitch_parts([OUT/f'{prefix}{split}_embedding_part_rank{r}.npz' for r in range(WORLD)],len(frame));np.save(cache,array)
   atomic_json(meta,{'data':identity,'rows':len(frame),'features':array.shape[1],'dtype':'float32','model':model_identity})
  barrier()
 del model
 if DEVICE.type=='cuda':torch.cuda.empty_cache()
 return True


class EarlyStopping:
 def __init__(self,patience):self.patience=patience;self.best=-float('inf');self.bad=0;self.stop=False
 def update(self,value):
  improved=value>self.best
  if improved:self.best=value;self.bad=0
  else:self.bad+=1
  self.stop=self.bad>=self.patience
  return improved


def loss_row(name,epoch,step,loss,validation=None):
 if RANK:return
 path=OUT/('smoke_losses.csv' if name.startswith('smoke_') else 'losses.csv');exists=path.exists()
 with path.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['model','epoch','global_step','training_loss','validation_mcc','time_utc'])
  if not exists:w.writeheader()
  w.writerow({'model':name,'epoch':epoch,'global_step':step,'training_loss':loss,'validation_mcc':validation,'time_utc':datetime.now(timezone.utc).isoformat()})


def baseline(config,deadline,smoke=False):
 if RANK:return
 prefix='smoke_' if smoke else '';seed_all(config['seed']);frames={s:load_frame(config[prefix+s+'_csv']) for s in ['train','validation']}
 arrays={};model_identity=source_identity(config)
 for split,frame in frames.items():
  meta=json.loads((OUT/f'{prefix}{split}_embeddings.json').read_text())
  if meta['data']!=data_identity(config[prefix+split+'_csv']):raise ValueError('Stale embeddings after split changed')
  require_model_identity(meta.get('model'),model_identity,f'{split} embedding cache')
  arrays[split]=np.load(OUT/f'{prefix}{split}_embeddings.npy',mmap_mode='r')
 head=torch.nn.Linear(arrays['train'].shape[1],1).to(DEVICE);optimizer=torch.optim.AdamW(head.parameters(),lr=config['baseline_learning_rate'],weight_decay=config['weight_decay'])
 stopper=EarlyStopping(config['early_stopping_patience']);step=0;best_path=OUT/f'{prefix}esm2_650m_frozen_best.pt';epochs=1 if smoke else config['baseline_epochs'];actual=0
 random_baseline(frames['validation'],config,prefix)
 for epoch in range(1,epochs+1):
  order=np.random.default_rng(config['seed']+epoch).permutation(len(frames['train']));total=0.;count=0;head.train()
  for start in range(0,len(order),config['baseline_batch_size']):
   if deadline_reached(deadline):break
   idx=order[start:start+config['baseline_batch_size']];x=torch.from_numpy(np.array(arrays['train'][idx])).to(DEVICE);y=torch.tensor(frames['train'].label.to_numpy()[idx],device=DEVICE,dtype=torch.float32)
   optimizer.zero_grad(set_to_none=True);loss=torch.nn.functional.binary_cross_entropy_with_logits(head(x).squeeze(-1),y)
   if not torch.isfinite(loss):raise FloatingPointError('Nonfinite baseline loss')
   loss.backward();torch.nn.utils.clip_grad_norm_(head.parameters(),config['grad_clip']);optimizer.step();step+=1;total+=loss.item()*len(idx);count+=len(idx)
  if not count:break
  head.eval();pred=[]
  with torch.inference_mode():
   for start in range(0,len(frames['validation']),8192):
    x=torch.from_numpy(np.array(arrays['validation'][start:start+8192])).to(DEVICE);pred.extend(torch.sigmoid(head(x).squeeze(-1)).cpu().tolist())
  m=metrics(frames['validation'].label,pred);actual=epoch;loss_row(prefix+'frozen',epoch,step,total/count,m['mcc']);log(f'baseline epoch={epoch} loss={total/count:.6f} validation={m}')
  if stopper.update(m['mcc']):
   state={'head.'+k:v.detach().cpu() for k,v in head.state_dict().items()}
   atomic_save(best_path,{'config':config,'state_dict':state,'epoch':epoch,'metrics':m,'source_checkpoint_identity':model_identity,'data_identity':{s:data_identity(config[prefix+s+'_csv']) for s in ['train','validation']},'training':{'mode':'frozen_embedding_linear_head','learning_rate':config['baseline_learning_rate'],'batch_size':config['baseline_batch_size'],'max_epochs':epochs,'seed':config['seed'],'warm_start':False}})
   write_predictions(prefix+'frozen',frames['validation'],pred,{'checkpoint':str(best_path),'best_epoch':epoch})
  if stopper.stop or deadline_reached(deadline):break
 if not best_path.exists():raise RuntimeError('No complete baseline validation/checkpoint was produced')
 atomic_json(OUT/f'{prefix}baseline_run.json',{'actual_epochs':actual,'global_steps':step,'early_stopped':stopper.stop,'deadline_reached':deadline_reached(deadline),'best_mcc':stopper.best,'head_architecture':'Linear(1280,1)','optimizer':'AdamW','learning_rate':config['baseline_learning_rate'],'max_epochs':epochs,'batch_size':config['baseline_batch_size']})


def evaluate_model(model,tokens,frame,config,deadline):
 model.eval();records=[];total=0.;count=0
 # Rotary caches survive evaluation and must remain usable by the next backward pass.
 with torch.no_grad():
  for idx in bucket_batches(tokens.lengths,config['evaluation_batch_size'],RANK,WORLD,config['seed'],False):
   if deadline_reached(deadline):break
   x,mask,y=tokens.batch(idx)
   with amp():logits=model(x,mask)
   if not torch.isfinite(logits).all():raise FloatingPointError('Nonfinite validation predictions')
   p=torch.sigmoid(logits.float()).cpu().tolist();records.extend(zip(idx,p))
 records=gather(records)
 if len(records)!=len(frame):return None
 records.sort();indices=[r[0] for r in records]
 if indices!=list(range(len(frame))):raise ValueError('Validation sharding duplicated or omitted rows')
 return np.asarray([r[1] for r in records])


def checkpoint(path,model,config,epoch,step,m,identity,extra):
 from modeling import trainable_state
 if RANK==0:atomic_save(path,{'config':config,'state_dict':trainable_state(model),'epoch':epoch,'global_step':step,'metrics':m,'source_checkpoint_identity':identity,**extra})


def restore_training_checkpoint(config,model,identity,expected_data,best_path,stopper,updates_per_epoch,prefix=''):
 path=Path(config['resume_checkpoint']);saved=torch.load(path,map_location='cpu',weights_only=False)
 require_model_identity(saved.get('source_checkpoint_identity'),identity,'training recovery')
 for key in ['model_dir','train_csv','validation_csv','seed','per_device_batch_size','expected_world_size','gradient_accumulation_steps','inject_layers','lora_r','lora_alpha','lora_dropout','gradient_checkpointing','learning_rate']:
  if saved['config'].get(key)!=config.get(key):raise ValueError(f'Recovery configuration changed: {key}')
 baseline_path=OUT/f'{prefix}esm2_650m_frozen_best.pt'
 if saved.get('training',{}).get('baseline_checkpoint_sha256')!=sha_file(baseline_path):raise ValueError('Recovery baseline/data provenance changed')
 if saved.get('data_identity',expected_data)!=expected_data:raise ValueError('Recovery data identity changed')
 step=int(saved['global_step'])
 if step<=0 or updates_per_epoch<=0:raise ValueError('Invalid recovery step or epoch size')
 completed,offset=divmod(step,updates_per_epoch)
 if int(saved['epoch'])!=(completed if offset==0 else completed+1):raise ValueError('Recovery epoch and global step disagree')
 best=torch.load(best_path,map_location='cpu',weights_only=False)
 require_model_identity(best.get('source_checkpoint_identity'),identity,'preserved best checkpoint')
 best_mcc=float(best['metrics']['mcc'])
 if not math.isfinite(best_mcc):raise ValueError('Invalid preserved best MCC')
 history={};history_path=OUT/(prefix+'losses.csv')
 if history_path.exists():
  with history_path.open() as f:
   for row in csv.DictReader(f):
    if row['model']==prefix+'lora' and row['validation_mcc']:
     epoch=int(row['epoch'])
     if int(row['global_step'])==epoch*updates_per_epoch and epoch<=completed:history[epoch]=float(row['validation_mcc'])
 if set(history)!=set(range(1,completed+1)):raise ValueError('Missing completed-epoch early-stopping history for recovery')
 for epoch in sorted(history):stopper.update(history[epoch])
 load_trainable_checkpoint(model,saved['state_dict'])
 recovery={'reason':'Rotary inference-cache failure repaired by no_grad validation','checkpoint':str(path),'checkpoint_sha256':sha_file(path),'from_global_step':step,'prior_epochs_completed':completed,'skipped_batches_in_resumed_epoch':offset*config.get('gradient_accumulation_steps',1),'optimizer_state_restored':False,'rng_state_restored':False,'note':'Original checkpoint contained adapter/head weights but no optimizer or RNG state. AdamW and RNG restart; deterministic batch order resumes after consumed batches. Resumed-epoch loss averages only post-recovery batches.'}
 return {'step':step,'start_epoch':completed+1,'completed_epochs':completed,'actual_epochs_started':int(saved['epoch']),'skip_batches':recovery['skipped_batches_in_resumed_epoch'],'best_mcc':best_mcc,'recovery':recovery}


def lora(config,deadline,smoke=False):
 from modeling import make_model
 if not smoke and WORLD!=config['expected_world_size']:raise ValueError('Full LoRA requires two DDP ranks; use torch.distributed.run --nproc_per_node=2')
 accum=config.get('gradient_accumulation_steps',1)
 if accum<1:raise ValueError('gradient_accumulation_steps must be positive')
 prefix='smoke_' if smoke else '';frames={s:load_frame(config[prefix+s+'_csv']) for s in ['train','validation']};tokens={s:Tokens(f,config['model_dir']) for s,f in frames.items()}
 seed_all(config['seed']);model=make_model(config,lora=True).to(DEVICE)
 baseline_path=OUT/f'{prefix}esm2_650m_frozen_best.pt'
 if not baseline_path.exists():raise FileNotFoundError(f'Fit frozen head first: {baseline_path}')
 frozen=torch.load(baseline_path,map_location='cpu',weights_only=False)
 require_model_identity(frozen.get('source_checkpoint_identity'),source_identity(config),'LoRA baseline warm start')
 expected_data={s:data_identity(config[prefix+s+'_csv']) for s in ['train','validation']}
 if frozen.get('data_identity')!=expected_data:raise ValueError('LoRA baseline warm start: stale or missing data identity; refit the baseline')
 model.head.load_state_dict({k.removeprefix('head.'):v for k,v in frozen['state_dict'].items() if k.startswith('head.')})
 identity=frozen['source_checkpoint_identity'];training={'warm_start':'Frozen-baseline trained linear head; LoRA adapters start from PEFT default initialization.','baseline_checkpoint':str(baseline_path),'baseline_checkpoint_sha256':sha_file(baseline_path),'optimizer':'AdamW','learning_rate':config['learning_rate'],'global_batch_size':config['per_device_batch_size']*WORLD*accum,'gradient_accumulation_steps':accum,'per_device_batch_size':config['per_device_batch_size'],'world_size':WORLD,'precision':'fp32 parameters, bf16 autocast on CUDA','drop_last':True,'initial_seed':config['seed']}
 if WORLD>1:
  from torch.nn.parallel import DistributedDataParallel
  wrapped=DistributedDataParallel(model,device_ids=[DEVICE.index] if DEVICE.type=='cuda' else None,broadcast_buffers=False,find_unused_parameters=False)
 else:wrapped=model
 params=[p for p in model.parameters() if p.requires_grad];optimizer=torch.optim.AdamW(params,lr=config['learning_rate'],weight_decay=config['weight_decay'])
 stopper=EarlyStopping(config['early_stopping_patience']);step=0;actual_epochs=0;completed_epochs=0;best_path=OUT/f'{prefix}esm2_650m_lora_best.pt';last_path=OUT/f'{prefix}esm2_650m_lora_latest.pt';epochs=1 if smoke else config['epochs'];stopped=False;best_pred=None;best_mcc=-float('inf');epoch_mcc=None;produced_best=False
 start_epoch=1;skip_batches=0
 if config.get('resume_checkpoint') and not smoke:
  updates_per_epoch=len(frames['train'])//(config['per_device_batch_size']*WORLD*accum)
  restored=restore_training_checkpoint(config,model,identity,expected_data,best_path,stopper,updates_per_epoch,prefix)
  step=restored['step'];start_epoch=restored['start_epoch'];skip_batches=restored['skip_batches'];completed_epochs=restored['completed_epochs'];actual_epochs=restored['actual_epochs_started'];best_mcc=restored['best_mcc'];produced_best=True
  training['recovery']=restored['recovery'];training['warm_start']='Continued LoRA and head weights from the archived checkpoint; original run started from the frozen baseline head.'
 training['data_identity']=expected_data
 log(f'LoRA trainable parameters={sum(p.numel() for p in params)} training={training}')
 random_baseline(frames['validation'],config,prefix)
 for epoch in range(start_epoch,epochs+1):
  batches=bucket_batches(tokens['train'].lengths,config['per_device_batch_size'],RANK,WORLD,config['seed']+epoch,True)
  batches=batches[:len(batches)//accum*accum]
  if not batches:raise ValueError('Training dataset smaller than global batch')
  full_batch_count=len(batches);batch_offset=skip_batches if epoch==start_epoch else 0
  batches=batches[batch_offset:]
  wrapped.train();epoch_loss=0.;epoch_steps=0;last_eval_step=-1
  for local_step,idx in enumerate(batches):
   if deadline_reached(deadline,collective=True):stopped=True;break
   x,mask,y=tokens['train'].batch(idx)
   if local_step%accum==0:optimizer.zero_grad(set_to_none=True)
   with amp():logits=wrapped(x,mask);loss=torch.nn.functional.binary_cross_entropy_with_logits(logits.float(),y)
   finite=torch.isfinite(loss).to(torch.int32)
   if WORLD>1:dist.all_reduce(finite,op=dist.ReduceOp.MIN)
   if not finite.item():raise FloatingPointError('Nonfinite training loss on at least one rank')
   (loss/accum).backward();epoch_steps+=1;epoch_loss+=loss.item()
   if (local_step+1)%accum:continue
   norm=torch.nn.utils.clip_grad_norm_(params,config['grad_clip'],error_if_nonfinite=True);optimizer.step();step+=1
   if step%config['log_every_steps']==0:
    log(f'lora epoch={epoch}/{epochs} batch={batch_offset+local_step+1}/{full_batch_count} step={step} loss={loss.item():.6f} grad_norm={float(norm):.4f}')
    logged_loss=loss.detach().float()
    if WORLD>1:dist.all_reduce(logged_loss);logged_loss/=WORLD
    loss_row(prefix+'lora',epoch,step,float(logged_loss))
   if step%config['checkpoint_every_steps']==0:checkpoint(last_path,model,config,epoch,step,None,identity,{'training':training,'complete_epoch':False})
   if step%config['evaluation_every_steps']==0 and local_step+1<len(batches):
    p=evaluate_model(model,tokens['validation'],frames['validation'],config,deadline);last_eval_step=step
    if p is None:stopped=True;break
    m=metrics(frames['validation'].label,p);log(f'interval validation step={step} metrics={m}');loss_row(prefix+'lora',epoch,step,epoch_loss/epoch_steps,m['mcc'])
    # Intermediate selection is allowed; patience counts only completed epochs.
    if m['mcc']>best_mcc:
     best_mcc=m['mcc'];best_pred=p;produced_best=True
     checkpoint(best_path,model,config,epoch,step,m,identity,{'training':training,'complete_epoch':False})
     write_predictions(prefix+'lora',frames['validation'],p,{'checkpoint':str(best_path),'best_epoch':epoch,'global_step':step,'selection':'maximum validation MCC'})
    wrapped.train()
  actual_epochs=epoch if epoch_steps else actual_epochs
  loss_total=torch.tensor([epoch_loss,epoch_steps],device=DEVICE,dtype=torch.float64)
  if WORLD>1:dist.all_reduce(loss_total)
  mean_loss=float(loss_total[0]/loss_total[1]) if loss_total[1]>0 else None
  complete=epoch_steps==len(batches)
  if complete:completed_epochs+=1
  checkpoint(last_path,model,config,epoch,step,None,identity,{'training':training,'complete_epoch':complete})
  if not stopped:
   p=evaluate_model(model,tokens['validation'],frames['validation'],config,deadline)
   if p is None:stopped=True
   else:
    m=metrics(frames['validation'].label,p);loss_row(prefix+'lora',epoch,step,mean_loss,m['mcc']);log(f'lora epoch={epoch} completed={complete} loss={mean_loss} validation={m}')
    stopper.update(m['mcc'])
    if m['mcc']>best_mcc:
     best_mcc=m['mcc'];best_pred=p;produced_best=True;checkpoint(best_path,model,config,epoch,step,m,identity,{'training':training,'complete_epoch':complete})
     write_predictions(prefix+'lora',frames['validation'],p,{'checkpoint':str(best_path),'best_epoch':epoch,'global_step':step,'selection':'maximum validation MCC'})
  if stopped or stopper.stop:break
 status='complete' if produced_best else ('incomplete_needs_final_validation' if step>0 else 'failed_no_training')
 if RANK==0:
  atomic_json(OUT/f'{prefix}lora_run.json',{'actual_epochs_started':actual_epochs,'completed_epochs':completed_epochs,'global_steps':step,'early_stopped':stopper.stop,'deadline_reached':stopped,'best_mcc':best_mcc if math.isfinite(best_mcc) else None,'best_checkpoint_exists':produced_best,'status':status,'latest_checkpoint':str(last_path),'selected_checkpoint':str(best_path) if produced_best else None,'training':training,'dropped_rows_per_epoch':len(frames['train'])%(config['per_device_batch_size']*WORLD*accum)})
 log(f'LoRA ended: completed_epochs={completed_epochs} steps={step} selected_checkpoint={produced_best}')
 if step==0:raise TimeoutError('LoRA deadline reached before any optimizer step; no training occurred')
 if not produced_best:log(f'Status incomplete_needs_final_validation: evaluate latest checkpoint {last_path} during the reserved final evaluation phase')
 return status


def load_trainable_checkpoint(model,state):
 expected={name for name,param in model.named_parameters() if param.requires_grad}
 missing_trainable=expected-set(state)
 if missing_trainable:raise ValueError(f'Checkpoint missing required trainable tensors: {sorted(missing_trainable)}')
 missing,unexpected=model.load_state_dict(state,strict=False)
 if unexpected:raise ValueError(f'Unexpected checkpoint weights: {unexpected}')


def standalone_evaluate(config,args,deadline):
 from modeling import make_model
 path=Path(args.checkpoint or OUT/'esm2_650m_lora_best.pt');saved=torch.load(path,map_location='cpu',weights_only=False)
 saved_config=dict(saved['config']);saved_config['model_dir']=config['model_dir']
 require_model_identity(source_identity(saved_config),saved.get('source_checkpoint_identity'),'saved checkpoint reload')
 use_lora=saved.get('training',{}).get('mode')!='frozen_embedding_linear_head'
 model=make_model(saved_config,lora=use_lora).to(DEVICE)
 load_trainable_checkpoint(model,saved['state_dict'])
 frame=load_frame(config[args.split+'_csv']);tokens=Tokens(frame,config['model_dir']);p=evaluate_model(model,tokens,frame,config,deadline)
 if p is None:raise TimeoutError('Evaluation deadline reached; partial metrics are not published')
 name=('lora' if use_lora else 'frozen')+('_test' if args.split=='test' else '')
 write_predictions(name,frame,p,{'checkpoint':str(path),'split':args.split,'checkpoint_sha256':sha_file(path),'data':data_identity(config[args.split+'_csv'])})


def resolve_config_paths(config, root=None):
    """Make path-valued config entries absolute against the project root.

    Paths are stored relative to the project root so that no absolute repository location is
    baked into the configuration.  Absolute entries are honoured unchanged.
    """
    base = Path(root).resolve() if root else ROOT
    for key in ("model_dir", "output_dir", "train_csv", "validation_csv", "test_csv",
                "smoke_train_csv", "smoke_validation_csv", "resume_checkpoint"):
        value = config.get(key)
        if isinstance(value, str) and value and not Path(value).is_absolute():
            config[key] = str((base / value).resolve())
    return config


def main():
 global OUT
 parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['extract','baseline','lora','smoke','evaluate'],required=True);parser.add_argument('--config',default=str(ROOT/'training_config.yaml'));parser.add_argument('--max-seconds',type=float,default=32400);parser.add_argument('--deadline',help='Absolute UTC ISO-8601 deadline');parser.add_argument('--checkpoint');parser.add_argument('--split',choices=['validation','test'],default='validation')
 args=parser.parse_args();config=resolve_config_paths(yaml.safe_load(Path(args.config).read_text()));OUT=Path(config['output_dir']);OUT.mkdir(parents=True,exist_ok=True)
 start=time.time();deadline=min(start+args.max_seconds,datetime.fromisoformat(args.deadline.replace('Z','+00:00')).timestamp()) if args.deadline else start+args.max_seconds
 setup();seed_all(config['seed']);log(f'Start mode={args.mode} world={WORLD} deadline={datetime.fromtimestamp(deadline,timezone.utc).isoformat()}')
 try:
  if args.mode=='extract':
   if not extraction(config,deadline):raise TimeoutError('Frozen embedding extraction incomplete at deadline')
  elif args.mode=='baseline':baseline(config,deadline)
  elif args.mode=='lora':lora(config,deadline)
  elif args.mode=='evaluate':standalone_evaluate(config,args,deadline)
  else:
   if extraction(config,deadline,True):
    baseline(config,deadline,True);barrier()
    if not deadline_reached(deadline,collective=True):lora(config,deadline,True)
    else:raise TimeoutError('Smoke deadline reached before LoRA stage')
   else:raise TimeoutError('Smoke embedding extraction incomplete')
  barrier();log(f'Finished mode={args.mode} elapsed_seconds={time.time()-start:.1f}')
 except BaseException:
  detail=traceback.format_exc();log(detail)
  with (OUT/f'error_rank{RANK}.log').open('a') as f:f.write(detail+'\n')
  with (OUT/'error_log.txt').open('a') as f:f.write(f'rank={RANK}\n{detail}\n')
  raise
 finally:
  if dist.is_initialized():dist.destroy_process_group()

if __name__=='__main__':main()
