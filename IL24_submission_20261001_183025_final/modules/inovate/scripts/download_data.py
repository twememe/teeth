import os,json,urllib.request,pathlib,hashlib,traceback,concurrent.futures,zipfile
ROOT=pathlib.Path(__file__).resolve().parents[1];os.chdir(ROOT);os.environ.update(json.loads((ROOT/'environment.json').read_text()))
RAW=ROOT/'data/raw';RAW.mkdir(parents=True,exist_ok=True)
def fetch(url,path):
 path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 if not path.exists():
  for attempt in range(2):
   try:
    with urllib.request.urlopen(url,timeout=180) as response:path.write_bytes(response.read())
    break
   except Exception:
    err=traceback.format_exc();print(err,flush=True)
    with (ROOT/'error_log.txt').open('a') as f:f.write(f'Download {url} attempt {attempt+1}\n'+err)
    if attempt:raise
 b=path.read_bytes();return {'url':url,'file':str(path.relative_to(ROOT)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def main():
 jobs=[];versions={}
 for short,repo in [('abbench','AbBibench/Antibody_Binding_Benchmark_Dataset'),('il6','alchemab/il6-binding-prediction')]:
  dest=RAW/short;dest.mkdir(exist_ok=True)
  fetch(f'https://huggingface.co/api/datasets/{repo}',dest/'repository.json')
  revision={'abbench':'556fd6913aa231c0d342a8658818be8f963fd582','il6':'f83d7fa0d08c143aefb4bb5b8a4e55ed0576c34f'}[short];versions[short]=revision
  fetch(f'https://huggingface.co/api/datasets/{repo}/tree/{revision}?recursive=true',dest/'tree.json')
  tree=json.loads((dest/'tree.json').read_text())
  for item in tree:
   if item['type']=='file' and (item['path'].endswith(('.csv','.parquet','.json','.md','.pdb'))):jobs.append((f'https://huggingface.co/datasets/{repo}/resolve/{revision}/'+item['path'],dest/item['path']))
 dest=RAW/'abagym';dest.mkdir(exist_ok=True)
 fetch('https://api.github.com/repos/3BioCompBio/AbAgym/git/trees/87af82ad68bb90921e15fb8b79c7bcb1d05b23d3?recursive=1',dest/'tree.json')
 tree=json.loads((dest/'tree.json').read_text());versions['abagym']=tree['sha']
 for item in tree['tree']:
  if item['path'] in ['AbAgym_data_non-redundant_interface.csv','AbAgym_metadata.csv','PDB_files.zip','README.md']:jobs.append((f"https://raw.githubusercontent.com/3BioCompBio/AbAgym/{tree['sha']}/"+item['path'],dest/item['path']))
 manifest=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for result in ex.map(lambda job:fetch(*job),jobs):manifest.append(result);print('Downloaded',result['file'],result['bytes'],flush=True)
 with zipfile.ZipFile(dest/'PDB_files.zip') as z:
  for member in z.infolist():
   target=(dest/member.filename).resolve()
   if not target.is_relative_to(dest.resolve()):raise ValueError('Unsafe archive path')
  z.extractall(dest)
 (ROOT/'data/download_manifest.json').write_text(json.dumps({'versions':versions,'files':manifest},indent=2))
 with (ROOT/'progress.txt').open('a') as f:f.write('All three requested datasets and official metadata/PDB downloaded with SHA256 manifest.\n')
if __name__=='__main__':main()
