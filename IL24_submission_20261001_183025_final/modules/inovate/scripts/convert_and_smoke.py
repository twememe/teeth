import os,pathlib,json,time,traceback,gc
os.chdir(pathlib.Path(__file__).resolve().parents[1])
os.environ.update(json.loads(pathlib.Path('environment.json').read_text()))
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import torch,esm
from transformers import EsmConfig,EsmModel,EsmTokenizer
from modeling import make_model,trainable_state

def main():
    torch.set_num_threads(8);torch.manual_seed(42)
    target=pathlib.Path('cache/esm2_650m_hf');target.mkdir(exist_ok=True)
    path=pathlib.Path('cache/torch/hub/checkpoints')
    # Trusted, previously acquired official Meta ESM checkpoint; never unpickle arbitrary data.
    native,alphabet=esm.pretrained.load_model_and_alphabet_core('esm2_t33_650M_UR50D',
        torch.load(path/'esm2_t33_650M_UR50D.pt',map_location='cpu',weights_only=False),
        torch.load(path/'esm2_t33_650M_UR50D-contact-regression.pt',map_location='cpu',weights_only=False))
    conf=EsmConfig(vocab_size=len(alphabet),hidden_size=1280,num_hidden_layers=33,num_attention_heads=20,
        intermediate_size=5120,pad_token_id=alphabet.padding_idx,mask_token_id=alphabet.mask_idx,
        max_position_embeddings=1026,position_embedding_type='rotary',emb_layer_norm_before=False,
        token_dropout=True,hidden_dropout_prob=0.,attention_probs_dropout_prob=0.,layer_norm_eps=1e-5)
    model=EsmModel(conf,add_pooling_layer=False);weights={}
    mapping={'embed_tokens.weight':'embeddings.word_embeddings.weight',
        'emb_layer_norm_after.weight':'encoder.emb_layer_norm_after.weight','emb_layer_norm_after.bias':'encoder.emb_layer_norm_after.bias',
        'contact_head.regression.weight':'contact_head.regression.weight','contact_head.regression.bias':'contact_head.regression.bias'}
    part={'self_attn.q_proj':'attention.self.query','self_attn.k_proj':'attention.self.key',
        'self_attn.v_proj':'attention.self.value','self_attn.out_proj':'attention.output.dense',
        'self_attn_layer_norm':'attention.LayerNorm','final_layer_norm':'LayerNorm',
        'fc1':'intermediate.dense','fc2':'output.dense','self_attn.rot_emb.inv_freq':'attention.self.rotary_embeddings.inv_freq'}
    expected=model.state_dict()
    for key,value in native.state_dict().items():
        dst=mapping.get(key)
        if key.startswith('layers.'):
            _,idx,sub=key.split('.',2)
            for old,new in part.items():
                if sub==old or sub.startswith(old+'.'):
                    dst='encoder.layer.'+idx+'.'+sub.replace(old,new,1);break
        if dst in expected:weights[dst]=value
    # HF allocates absolute positions even for rotary ESM; these are never used.
    weights['embeddings.position_embeddings.weight']=torch.zeros_like(expected['embeddings.position_embeddings.weight'])
    missing=set(expected)-set(weights)
    assert not missing,missing
    model.load_state_dict(weights,strict=True)
    (target/'vocab.txt').write_text('\n'.join(alphabet.all_toks)+'\n')
    tokenizer=EsmTokenizer(str(target/'vocab.txt'));tokenizer.save_pretrained(target)
    seqs=['QVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGYYAMDYWGQGTLVTVSS',
          'DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPLTFGQGTKVEIK']
    batch=tokenizer(seqs,return_tensors='pt',padding=True);device=torch.device('cuda:1')
    native=native.to(device).eval();model=model.to(device).eval();batch={k:v.to(device) for k,v in batch.items()}
    with torch.no_grad():
        a=native(batch['input_ids'],repr_layers=[33])['representations'][33]
        b=model(**batch).last_hidden_state
    valid=batch['attention_mask'].bool();diff=(a[valid]-b[valid]).abs().max().item()
    print('CONVERSION_MAX_ABS_ERROR',diff,flush=True)
    torch.testing.assert_close(a[valid],b[valid],atol=2e-4,rtol=2e-4)
    model.cpu().save_pretrained(target,safe_serialization=True)
    del native,model,a,b,weights,expected;gc.collect();torch.cuda.empty_cache()
    cfg={'model_dir':str(target.resolve()),'sdpa':True,'lora_r':8,'lora_alpha':16,'lora_dropout':0.1,
        'inject_layers':list(range(16)),'gradient_checkpointing':False}
    m=make_model(cfg,lora=True).to(device).train()
    trainable=[n for n,p in m.named_parameters() if p.requires_grad]
    assert all('lora_' in n or n.startswith('head.') for n in trainable)
    opt=torch.optim.AdamW((p for p in m.parameters() if p.requires_grad),lr=1e-4)
    y=torch.tensor([0.,1.],device=device);losses=[];start=time.time()
    first={n:p.detach().clone() for n,p in m.named_parameters() if p.requires_grad}
    for step in range(3):
        opt.zero_grad(set_to_none=True)
        with torch.autocast('cuda',dtype=torch.bfloat16):
            logits=m(**batch);loss=torch.nn.functional.binary_cross_entropy_with_logits(logits.float(),y)
        assert torch.isfinite(loss);loss.backward()
        assert any(p.grad is not None and torch.isfinite(p.grad).all() and p.grad.abs().max()>0 for n,p in m.named_parameters() if 'lora_B' in n)
        torch.nn.utils.clip_grad_norm_((p for p in m.parameters() if p.requires_grad),1.)
        opt.step();losses.append(loss.item());print('SMOKE_STEP',step,'LOSS',loss.item(),flush=True)
    changed=[n for n,p in m.named_parameters() if n in first and not torch.equal(first[n],p)]
    assert any('lora_B' in n for n in changed)
    checkpoint={'state_dict':trainable_state(m),'config':cfg,'kind':'smoke_only_not_final','steps':3}
    torch.save(checkpoint,'output/esm2_650m_smoke.pt')
    m.eval()
    with torch.no_grad(),torch.autocast('cuda',dtype=torch.bfloat16):pred=m(**batch).float().cpu()
    sd=torch.load('output/esm2_650m_smoke.pt',map_location='cpu',weights_only=True)['state_dict']
    for n,p in m.named_parameters():
        if p.requires_grad:p.data.zero_()
    m.load_state_dict(sd,strict=False)
    with torch.no_grad(),torch.autocast('cuda',dtype=torch.bfloat16):restored=m(**batch).float().cpu()
    torch.testing.assert_close(pred,restored,rtol=0,atol=0)
    report={'status':'passed','backbone':'ESM-2 650M','backend':'PEFT + Transformers','gpu':'NVIDIA RTX A6000 GPU1',
        'real_forward_backward_steps':3,'losses':losses,'seconds':time.time()-start,'trainable_parameters':sum(p.numel() for p in m.parameters() if p.requires_grad),
        'max_gpu_allocated_gib':torch.cuda.max_memory_allocated(device)/1024**3,'native_conversion_max_abs_error':diff,
        'adapter_gradients_verified':True,'checkpoint_reload_verified':True,'not_performance_evaluation':True}
    pathlib.Path('output/gpu_smoke_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
    with open('progress.txt','a') as f:f.write(time.strftime('%Y-%m-%d %H:%M:%S ')+'ESM-2 650M real GPU forward/backward/save/reload feasibility passed.\n')

if __name__=='__main__':
    try:main()
    except Exception:
        err=traceback.format_exc();print(err,flush=True)
        with open('error_log.txt','a') as f:f.write(err)
        raise
