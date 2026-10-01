# 冻结基座外置与等价验证

本项目实际为“官方ESM-2基座冻结、本项目训练LoRA适配器/线性分类头”。最终checkpoint含64个LoRA张量和2个分类头张量，没有合并基座；全部3份原有项目checkpoint保留，不按扩展名批量排除权重。

官方原始格式是fair-esm .pt，正式分类器使用Transformers EsmModel，故沿用原转换脚本中的必要部分：参数改名映射、加载contact-regression、保留float32、生成原有零值绝对位置表（rotary不使用）、不装入原分类器不使用的官方masked-LM head。原脚本后面的训练smoke步骤不执行。不量化、不重训、不改变模型和筛选规则。

## 固定来源

官方仓库 https://github.com/facebookresearch/esm 。加载器fair-esm 2.0.0（v2.0.0）完整提交为`0b59d87ebef95948c735b1f7aad463dc6dfa991b`。权重文件通过官方Meta对象存储发行，完整对象版本如下；这与加载器Git提交是两个独立版本。历史未记录权重Git revision，不伪造历史revision。本次发布固定并验证能重建A的发行文件内容，URL携带versionId且下载必须同时满足完整SHA。

| 官方文件 | 完整S3对象版本 | 文件SHA256 |
|---|---|---|
| esm2_t33_650M_UR50D.pt | 2G4typsUSwOKQLH35sqPMJ27ZLX1AEJG | ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c |
| esm2_t33_650M_UR50D-contact-regression.pt | cEcpCgy3_xlt3NCWsoXOzeIrh7.EDOBu | 8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575 |

准确下载URL、文件字节数、转换依赖和实际配置见`models/base_manifest.json`。官方文件只与自己的源SHA比较。转换文件的SHA不同于官方.pt是正常现象，不能据此判模型错误。

## 独立参考与指纹

先从未修改的原正式基座加载参考A，再从原始官方文件在独立目录准备B；不是B保存再重载后自证。固定A记录位于`models/base_reference.json`，准备程序没有生成或覆盖它的功能。固定参考文件自身也受manifest与发布SHA清单约束。

按参数名固定升序，覆盖state_dict及全部注册buffer（包括非持久position_ids）。每项包含完整名字、torch dtype、shape及原dtype连续CPU字节；小端内容，无任何降精度。SHA起始前缀为`IL24-ESM2-TENSORS-v1`加NUL；每项依次加入8字节大端header长度、规范JSON header、8字节大端内容长度、张量内容。JSON按键排序、无多余空格、ensure_ascii=True。别名按每个名字重复纳入并单独比较共享组；实际A无共享组。配置仅排除位置字段_name_or_path；词表、特殊token/ID和tokenizer行为也独立校验。实际加载严格检查缺失、额外、shape冲突键，AutoTokenizer必须选用EsmTokenizer。

旧项目checkpoint仍原字节保留，其历史config/文件SHA身份只映射到这个已验证的A参考；新运行先验证实际基座张量，再比较身份。不会将B的文件SHA伪装成A的文件SHA，也不会接受任意未知旧身份。

## 上一轮外置基座等价结果

568个张量的键、shape、dtype、值全部精确相同，最大张量差0；加载配置、tokenizer及8样例实际token序列相同。使用同一份项目checkpoint、同一A6000、同一软件、原eval/no_grad/SDPA/bf16 autocast/分桶/seed42，8条真实样例的原始logit最大差0、概率最大差0、阈值0.5标签全部相同。两次正式run.py输出与上一版参考CSV逐字节一致；容差rtol=atol=0。一次受控张量改动被拒绝，且不改磁盘上的模型。

记录：`logs/base_externalization/equivalence.json`、`A_B_scores.csv`、`inference_A.log`、`inference_B.log`。转换文件字节SHA在本环境也恰好一致，但不是验收条件；源文件SHA与转换文件SHA不同，各自意义已分开。后续使用者只需固定清单与随包8样例，不需要A的大文件。

## 实际准备与推理

完整可复制命令见README的“冻结基座准备”。`prepare_base.py`无参数从固定官方URL获取且复用已校验缓存；`--official-dir ./official_files`只读离线原始.pt；`--converted-dir ./offline_model --verify-only`只读校验转换模型，再用`run.py --model-dir ./offline_model`加载。原始与转换输入类型不能混用。转换写入输出父目录下的独立临时目录，完成校验后改名，不修改官方缓存或原工程。

只下载两个必要官方文件，共2,604,541,236字节；运行模型约2.603 GB，额外磁盘建议6 GB、CPU内存建议16 GB。转换用CPU；推理推荐bf16 GPU。原软件版本锁定于manifest/esm2.txt，安装顺序见envs/README。历史外置基座验证复用现有环境，未重新安装环境或训练。自动路径实际测试复用了大文件缓存、下载了缺少的官方回归文件，并成功验证模型；不重复下载大权重。

ZIP只移出明确识别的冻结运行基座目录；项目训练权重、原配置/tokenizer和参考记录均保留。外置基座/缓存不属于静态SHA清单，由准备入口校验；新ZIP的解压验证结果及ZIP摘要写同级sidecar。

本次最终收尾的实际新环境安装、空缓存官方下载、CPU转换及GPU推理结论只见ZIP同级 `IL24_submission_20261001_183025_final.validation.json` / `.validation.md`，不由以上历史记录代替。首次执行顺序统一为环境检查 → `python tools/prepare_base.py` → `python run.py --module esm2 --output ./base_model_check`。
