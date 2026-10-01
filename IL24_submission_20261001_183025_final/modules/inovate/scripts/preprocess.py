"""Versioned-source preprocessing. Plain H+L concatenation and observed PDB-chain mutants."""
import csv,json,hashlib,math,random,re,statistics,traceback,collections,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data';AA=set('ACDEFGHIKLMNPQRSTVWY');SEED=42
AA3=dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),'ACDEFGHIKLMNPQRSTVWY'))
def valid_label(value):
 try:v=float(value)
 except (ValueError,TypeError):raise ValueError('invalid_label')
 if not math.isfinite(v) or v not in (0,1):raise ValueError('invalid_label')
 return int(v)
def valid_sequence(seq):
 if not seq or set(seq)-AA:raise ValueError('invalid_amino_acids')
 if len(seq)>1022:raise ValueError('sequence_length_above_1022')
 return seq

def mutate_sequence(wt,notation,residue_map=None):
 valid_sequence(wt);chars=list(wt);used=set()
 for token in notation.split(':'):
  match=re.fullmatch(r'([ACDEFGHIKLMNPQRSTVWY])(\d+)([A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])',token)
  if not match:raise ValueError('invalid_mutation_notation')
  before,number,icode,after=match.groups()
  if residue_map is None:
   if icode:raise ValueError('insertion_code_requires_pdb_mapping')
   index=int(number)-1
  else:
   if (number,icode) not in residue_map:raise ValueError('mutation_site_missing_in_pdb')
   index=residue_map[(number,icode)]
  if index<0 or index>=len(wt):raise ValueError('mutation_out_of_bounds')
  if wt[index]!=before:raise ValueError('wildtype_mismatch')
  if index in used:raise ValueError('duplicate_mutation_position')
  used.add(index);chars[index]=after
 return ''.join(chars)

def pdb_chains(path):
 chains={};seen={}
 for line in path.read_text().splitlines():
  if line.startswith('ENDMDL'):break
  if not line.startswith('ATOM  ') or line[12:16].strip()!='CA' or line[16] not in (' ','A'):continue
  chain=line[21];key=(line[22:26].strip(),line[26].strip());res=AA3.get(line[17:20].strip(),'X')
  if key in seen.setdefault(chain,set()):continue
  seen[chain].add(key);chains.setdefault(chain,[]).append((key,res))
 return {chain:(''.join(x[1] for x in residues),{x[0]:i for i,x in enumerate(residues)}) for chain,residues in chains.items()}

def split_group_map(rows,train_target=200000,val_target=20000):
 groups=collections.defaultdict(list)
 for row in rows:groups[row['model_sequence']].append(row)
 official={};priority={'train':0,'validation':1,'test':2}
 for seq,rs in groups.items():
  splits={r.get('official_split','') for r in rs}-{''}
  if splits:official[seq]=max(splits,key=lambda s:priority[s])
 n=len(rows)
 if n<train_target+val_target+max(1,int(.01*n)):
  train_target=int(n*.8);val_target=int(n*.1)
 assignment=dict(official);counts=collections.Counter()
 for seq,split in official.items():counts[split]+=len(groups[seq])
 keys=sorted(set(groups)-set(official));random.Random(SEED).shuffle(keys)
 for seq in keys:
  count=len(groups[seq])
  split='train' if counts['train']+count<=train_target else ('validation' if counts['validation']+count<=val_target else 'test')
  assignment[seq]=split;counts[split]+=count
 return assignment

def validate_splits(splits):
 seen=set()
 for name,rows in splits.items():
  seqs={r['model_sequence'] for r in rows}
  assert not seen.intersection(seqs),f'sequence leakage: {name}'
  seen.update(seqs)
  for r in rows:
   valid_sequence(r['model_sequence'])
   if 'label' in r:valid_label(r['label'])
 return True

def write_csv(path,rows,fields=None):
 if fields is None:fields=list(dict.fromkeys(k for row in rows for k in row))
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def main():
 import pandas as pd
 candidates=[];rejects=[];raw_counts={};ab_rows=[];ag_rows=[];il_rows=[]
 def reject(row,reason):rejects.append({'row_id':row.get('row_id',''),'source':row.get('source',''),'reason':reason,'model_sequence':row.get('model_sequence','')})
 def add(row):
  try:
   valid_sequence(row['model_sequence'])
   if not math.isfinite(float(row['raw_score'])):raise ValueError('nonfinite_score')
   if row['source']!='abbench':row['label']=valid_label(row['label'])
   row['group_id']=hashlib.sha256(row['model_sequence'].encode()).hexdigest();row['record_id']=row['row_id'];row['original_score']=row['raw_score'];candidates.append(row)
  except (ValueError,TypeError) as e:reject(row,str(e))
 for p in sorted((DATA/'raw/abbench/binding_affinity').glob('*.csv')):
  for i,x in enumerate(csv.DictReader(p.open())):
   h=x.get('heavy_chain_seq','').strip().upper();l=x.get('light_chain_seq','').strip().upper();rid=f'abbench:{p.stem}:{i}'
   out=dict(x,heavy_sequence=h,light_sequence=l,mutation='',affinity_score=x['binding_score'],source_file=p.name,row_id=rid,mutation_status='not_provided_by_source')
   ab_rows.append(out)
   add(dict(sequence=h+l,model_sequence=h+l,mutation='',source='abbench',raw_score=x['binding_score'],row_id=rid,official_split='',assay=p.stem,mutated_partner='antibody',label_semantics='above_training_median_numeric_binding_score',representation='heavy_then_light_concatenation_without_linker',mutation_status='not_provided_by_source',heavy_sequence=h,light_sequence=l))
 write_csv(DATA/'abbench_full.csv',ab_rows);raw_counts['abbench']=len(ab_rows)
 metadata={x['DMS_name']:x for x in csv.DictReader((DATA/'raw/abagym/AbAgym_metadata.csv').open())}
 pdbs={p.stem:p for p in (DATA/'raw/abagym').rglob('*.pdb')};cache={}
 for i,x in enumerate(csv.DictReader((DATA/'raw/abagym/AbAgym_data_non-redundant_interface.csv').open())):
  meta=metadata[x['DMS_name']];rid=f'abagym:{i}';score=float(x['MinMax_normalized_DMS_score']);ag_rows.append(dict(x,DMS_score_MinMax=x['MinMax_normalized_DMS_score'],row_id=rid,**{k:v for k,v in meta.items() if k not in x}))
  row=dict(sequence='',model_sequence='',mutation=x['mut_names'],source='abagym',raw_score=score,row_id=rid,official_split='',assay=x['DMS_name'],mutated_partner=meta['DMS_on'],label_semantics='high_numeric_MinMax_DMS_score',experimental_DMS_type=meta['experimental_DMS_type'],PDB_file=x['PDB_file'],wildtype=x['wildtype'],site=x['site'],chains=x['chains'],DMS_score=x['DMS_score'],representation='representative_mutated_chain_PDB_observed_residues')
  try:
   if score==.5:raise ValueError('threshold_tie_0.5')
   if not 0<=score<=1:raise ValueError('MinMax_out_of_range')
   if x['PDB_file'] not in pdbs:raise ValueError('pdb_not_found')
   if x['PDB_file'] not in cache:cache[x['PDB_file']]=pdb_chains(pdbs[x['PDB_file']])
   structure=cache[x['PDB_file']];verified=[];notation=x['wildtype']+x['site']+x['mutation']
   for chain in x['chains']:
    if chain not in structure:raise ValueError('chain_not_found')
    wt,mapping=structure[chain];mutant=mutate_sequence(wt,notation,mapping);verified.append((len(wt),chain,wt,mutant))
   _,chain,wt,mutant=max(verified,key=lambda t:(t[0],t[1]))
   row.update(sequence=wt+':'+notation,model_sequence=mutant,wildtype_sequence=wt,representative_chain=chain,mutation_notation=notation,label=int(score>.5),mutation_status='verified_PDB_numbering')
   add(row)
  except ValueError as e:reject(row,str(e))
 write_csv(DATA/'abagym_interface.csv',ag_rows);raw_counts['abagym']=len(ag_rows)
 for p in sorted((DATA/'raw/il6/data').glob('*.parquet')):
  split={'eval':'validation','train':'train','test':'test'}[p.name.split('-')[0]]
  for i,x in enumerate(pd.read_parquet(p).to_dict('records')):
   rid=f'il6:{split}:{i}';seq=x['vh'].strip().upper();il_rows.append(dict(x,sequence=seq,row_id=rid,official_split=split))
   add(dict(sequence=seq,model_sequence=seq,mutation='',label=x['label'],source='il6',raw_score=x['label'],row_id=rid,official_split=split,assay='IL6_binding',mutated_partner='antibody',label_semantics='original_nonbinder_0_binder_1',representation='source_vh',mutation_status='not_provided_by_source'))
 write_csv(DATA/'il6.csv',il_rows);raw_counts['il6']=len(il_rows)
 assignments=split_group_map(candidates);train_scores=[float(r['raw_score']) for r in candidates if r['source']=='abbench' and assignments[r['model_sequence']]=='train'];median=statistics.median(train_scores)
 splits={s:[] for s in ['train','validation','test']}
 for row in candidates:
  if row['source']=='abbench':
   score=float(row['raw_score'])
   if score==median:reject(row,'training_median_tie');continue
   row['label']=int(score>median)
  splits[assignments[row['model_sequence']]].append(row)
 validate_splits(splits)
 fields=['sequence','model_sequence','mutation','label','source','raw_score','group_id','row_id','record_id','original_score','official_split','assay','mutated_partner','label_semantics','representation','mutation_status','heavy_sequence','light_sequence','experimental_DMS_type','PDB_file','wildtype','site','chains','DMS_score','wildtype_sequence','representative_chain','mutation_notation']
 for split,rows in splits.items():
  random.Random(SEED).shuffle(rows);write_csv(DATA/f'{split}.csv',rows,fields)
 for split,n in [('train',1000),('validation',200)]:
  # Representative deterministic smoke subset, including all three sources where available.
  bysource={s:[r for r in splits[split] if r['source']==s] for s in raw_counts};chosen=[]
  quota=min(50,n//len(bysource))
  for rs in bysource.values():chosen.extend(rs[:quota])
  taken={r['row_id'] for r in chosen};chosen.extend(r for r in splits[split] if r['row_id'] not in taken);chosen=chosen[:n]
  random.Random(SEED).shuffle(chosen);write_csv(DATA/f'smoke_{split}.csv',chosen,fields)
 write_csv(DATA/'rejected_records.csv',rejects,['row_id','source','reason','model_sequence'])
 allrows=sum(splits.values(),[]);group_labels=collections.defaultdict(set)
 for row in allrows:group_labels[row['group_id']].add(row['label'])
 audit={'seed':SEED,'python':platform.python_version(),'pandas':pd.__version__,'raw_counts':raw_counts,'total_raw':sum(raw_counts.values()),'candidate_counts':dict(collections.Counter(r['source'] for r in candidates)),'rejections':dict(collections.Counter(r['source']+':'+r['reason'] for r in rejects)),'abbench_training_median':median,'abbench_training_median_fit_rows':len(train_scores),'threshold_fit_before_tie_exclusion':True,'requested_counts':{'train':200000,'validation':20000},'split_counts':{s:len(rs) for s,rs in splits.items()},'unique_model_sequences':len(group_labels),'conflicting_label_sequence_groups':sum(len(v)>1 for v in group_labels.values()),'split_unique_sequences':{s:len({r['model_sequence'] for r in rs}) for s,rs in splits.items()},'source_split_label_counts':{s:dict(collections.Counter(f"{r['source']}|{r['label']}" for r in rs)) for s,rs in splits.items()},'partner_split_counts':{s:dict(collections.Counter(f"{r['source']}|{r['mutated_partner']}" for r in rs)) for s,rs in splits.items()},'official_il6_split_changes':sum(r['source']=='il6' and r['official_split']!=s for s,rs in splits.items() for r in rs),'max_model_length':max(map(lambda r:len(r['model_sequence']),allrows)),'checks':{'identical_sequence_splits_disjoint':True,'canonical_AA_only':True,'labels_binary':True,'all_model_lengths_at_most_1022':True},'split_sha256':{s:hashlib.sha256((DATA/f'{s}.csv').read_bytes()).hexdigest() for s in splits},'manifest':'data/download_manifest.json','limitations':['Pooled labels are heterogeneous; high numeric DMS score can represent immune escape, not antibody affinity improvement.','AbBiBench mutation annotations are absent from source; mutation stays empty. Source sequences are already measured variants.','AbBiBench H and L are concatenated without linker; this introduces an artificial chain boundary.','AbAgym model input is the longest verified representative mutated PDB chain with observed ATOM residues only, not a complete antibody-antigen complex. Unresolved residues are absent. All equivalent chains mutation sites must verify.','Exact duplicate sequences are grouped across all sources; sequence identity/antigen/family-disjoint generalization is not assessed. Conflicting assays/labels are retained for feasibility and reported.','Original IL6 split assignments are retained, with test-over-validation-over-train priority if an identical sequence has conflicting official assignments.','No oversampling or fabricated records. Split targets apply to eligible rows before AbBiBench median ties; deficits are reported.']}
 (ROOT/'output/data_audit.json').write_text(json.dumps(audit,indent=2))
 report='# Data implementation report\n\n'+json.dumps({k:audit[k] for k in ['raw_counts','total_raw','split_counts','unique_model_sequences','rejections','abbench_training_median','official_il6_split_changes']},indent=2)+'\n\n'+'\n'.join('- '+s for s in audit['limitations'])+'\n\nSource revisions and hashes: data/download_manifest.json. Rejected rows: data/rejected_records.csv. Smoke files contain 1,000 training and 200 validation rows with all available sources represented.\n'
 (ROOT/'output/data_implementation_report.md').write_text(report)
 print(json.dumps(audit,indent=2),flush=True)
 with (ROOT/'progress.txt').open('a') as f:f.write(f'Data preprocessing complete: {audit["split_counts"]}; raw {audit["raw_counts"]}; median {median}; disjoint canonical model sequences verified.\n')
if __name__=='__main__':
 try:main()
 except Exception:
  error=traceback.format_exc();print(error,flush=True)
  with (ROOT/'error_log.txt').open('a') as f:f.write(error)
  raise
