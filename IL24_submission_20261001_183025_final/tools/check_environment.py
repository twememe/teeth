#!/usr/bin/env python3
import argparse,sys,platform,json,importlib,importlib.metadata
ap=argparse.ArgumentParser();ap.add_argument('--module',choices=['export','postprocess','phase2','esm2'],required=True);a=ap.parse_args()
names={'export':[],'postprocess':['numpy','pandas','scipy','matplotlib','Bio','yaml'],'phase2':['numpy','pandas','scipy','Bio','torch','ImmuneBuilder','prodigy_prot','boltz'],'esm2':['torch','transformers','peft','safetensors','numpy','pandas','scipy','sklearn','yaml']}
result={'module':a.module,'python':sys.version,'platform':platform.platform(),'imports':{}}
for name in names[a.module]:
 m=importlib.import_module(name);result['imports'][name]={'version':getattr(m,'__version__','see distribution metadata'),'file':getattr(m,'__file__',None)}
if 'torch' in names[a.module]:
 import torch
 result['torch_cuda']=torch.version.cuda;result['cuda_available']=torch.cuda.is_available();result['gpu']=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
print(json.dumps(result,indent=2));print('ENVIRONMENT_IMPORT_CHECK_PASSED')
