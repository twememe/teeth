"""ESM-2 sequence classifier: frozen encoder or PEFT LoRA with the same linear head."""
import math, types
import torch
from torch import nn
from torch.nn import functional as F
from transformers import EsmModel
from peft import LoraConfig, get_peft_model

def add_lora(backbone,layers,r=8,alpha=16,dropout=0.1):
    return get_peft_model(backbone,LoraConfig(r=r,lora_alpha=alpha,lora_dropout=dropout,
        target_modules=['query','value'],layers_to_transform=list(layers),layers_pattern='layer',bias='none'))

def _sdpa_forward(self,hidden_states,attention_mask=None,head_mask=None,encoder_hidden_states=None,
                  encoder_attention_mask=None,past_key_value=None,output_attentions=False):
    if head_mask is not None or encoder_hidden_states is not None or past_key_value is not None or output_attentions or self.is_decoder:
        return self._original_forward(hidden_states,attention_mask,head_mask,encoder_hidden_states,encoder_attention_mask,past_key_value,output_attentions)
    q=self.transpose_for_scores(self.query(hidden_states))*self.attention_head_size**-0.5
    k=self.transpose_for_scores(self.key(hidden_states));v=self.transpose_for_scores(self.value(hidden_states))
    if self.position_embedding_type=='rotary':q,k=self.rotary_embeddings(q,k)
    out=F.scaled_dot_product_attention(q,k,v,attn_mask=attention_mask,
        dropout_p=self.dropout.p if self.training else 0.,scale=1.0)
    out=out.permute(0,2,1,3).contiguous().view(hidden_states.shape[0],hidden_states.shape[1],self.all_head_size)
    return (out,)

def enable_sdpa(backbone):
    from transformers.models.esm.modeling_esm import EsmSelfAttention
    for m in backbone.modules():
        if isinstance(m,EsmSelfAttention) and not hasattr(m,'_original_forward'):
            m._original_forward=m.forward
            m.forward=types.MethodType(_sdpa_forward,m)

class SequenceClassifier(nn.Module):
    def __init__(self,backbone):
        super().__init__();self.backbone=backbone
        self.head=nn.Linear(backbone.config.hidden_size,1)
    def features(self,input_ids,attention_mask):
        hidden=self.backbone(input_ids=input_ids,attention_mask=attention_mask,return_dict=True).last_hidden_state
        valid=(attention_mask.bool() & input_ids.ne(0) & input_ids.ne(1) & input_ids.ne(2)).unsqueeze(-1)
        return (hidden.float()*valid).sum(1)/valid.sum(1).clamp_min(1)
    def forward(self,input_ids,attention_mask,features_only=False):
        features=self.features(input_ids,attention_mask)
        return features if features_only else self.head(features).squeeze(-1)

def make_model(config,lora=True):
    backbone=EsmModel.from_pretrained(config['model_dir'],add_pooling_layer=False,local_files_only=True)
    for p in backbone.parameters():p.requires_grad=False
    if config.get('sdpa',True):enable_sdpa(backbone)
    if lora:
        backbone=add_lora(backbone,config.get('inject_layers',list(range(16))),config.get('lora_r',8),config.get('lora_alpha',16),config.get('lora_dropout',0.1))
        if config.get('gradient_checkpointing',False):
            backbone.gradient_checkpointing_enable(gradient_checkpointing_kwargs={'use_reentrant':False})
    return SequenceClassifier(backbone)

def trainable_state(model):
    names={n for n,p in model.named_parameters() if p.requires_grad}
    return {k:v.detach().cpu() for k,v in model.state_dict().items() if k in names}
