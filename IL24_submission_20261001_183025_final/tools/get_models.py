#!/usr/bin/env python3
"""Acquire recorded official assets and require the project's frozen SHA256 identity."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--group',choices=['boltz','antibody','phase1'],required=True);ap.add_argument('--destination',type=Path,required=True);ap.add_argument('--verify-only',action='store_true');a=ap.parse_args()
 manifest='phase1_assets.json' if a.group=='phase1' else 'external_assets.json'
 assets=[x for x in json.loads((ROOT/'models'/manifest).read_text()) if x['group']==a.group]
 for asset in assets:
  p=a.destination/asset['asset'];p.parent.mkdir(parents=True,exist_ok=True)
  if not p.exists():
   if a.verify_only:raise FileNotFoundError(p)
   temp=p.with_suffix(p.suffix+'.partial')
   with urllib.request.urlopen(asset['url'],timeout=180) as src,temp.open('wb') as out:
    for block in iter(lambda:src.read(8*1024*1024),b''):out.write(block)
   if temp.stat().st_size!=asset['bytes'] or sha(temp)!=asset['sha256']:raise ValueError('Downloaded identity mismatch: '+str(p))
   temp.rename(p)
  if p.stat().st_size!=asset['bytes'] or sha(p)!=asset['sha256']:raise ValueError('Asset identity mismatch: '+str(p))
  print('SHA256 OK',asset['asset'])
if __name__=='__main__':main()
