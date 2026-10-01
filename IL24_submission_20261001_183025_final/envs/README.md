# 环境安装

默认导出只需Python>=3.10标准库。CPU轻量复核首选：

```bash
python3.10 -m venv .env-post
. .env-post/bin/activate
python -m pip install -r envs/postprocess.txt
python -m pip check
python tools/check_environment.py --module postprocess
python run.py --module phase1 --output ./phase1_check
python run.py --module phase2 --output ./phase2_check
```

ESM首选Linux x86_64、Python3.10，NVIDIA GPU支持bf16，CUDA12.4构建PyTorch；历史环境驱动570.133.07、A6000 48GB。安装顺序：

```bash
python3.10 -m venv .env-esm2
. .env-esm2/bin/activate
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
python -m pip install -r envs/esm2.txt
python -m pip check
python tools/check_environment.py --module esm2
python tools/prepare_base.py
python run.py --module esm2 --output ./esm2_check
```

完整Phase2（同Python/CUDA基础，单独环境）：

```bash
python3.10 -m venv .env-phase2
. .env-phase2/bin/activate
python -m pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 --index-url https://download.pytorch.org/whl/cu124
python -m pip install -r envs/phase2.txt
python -m pip install --no-deps ./modules/phase2/vendor/boltz
python tools/check_environment.py --module phase2
```

`torch==2.6.0+cu124`是实测版本标识；不能先从默认源安装另一torch后假称相同构建。PRODIGY实际分发名为`prodigy-prot`，不是误写的PRODIGY包。HMMER3.4随包提供原工具；ANARCI、OpenMM的原依赖在Phase2清单中。Full181源代码及MSA runtime代码随包，不使用原服务器editable路径。Linux glibc2.35/kernel6.8为历史记录；其他平台未验证。

Phase1上游分别用原`modules/phase1ab/envs`锁定规格与`envs/phase1_discotope.txt`：核心含MAFFT/FreeSASA2.2.1/DSSP4.6.1（libmcfp2.0.1）；BepiPred原torch2.6+cu124、bp3 0.0.12.7/fair-esm1.0.3；DiscoTope原Python3.14.7、torch2.10.0 CUDA12.8。不能强行将这些都安装进ESM环境。历史核心平台锁定文件只是原版本证据，前缀已从可用YAML中移除。

历史发布复用了旧环境。2026-10-01最终收尾的新CPU/ESM环境安装命令、版本、退出码与结论在ZIP同级验收sidecar中，以实际记录为准；完整上游研究工具不在本轮重装范围。历史精确实测包见`logs/validation/runtime_*.json`。CUDA Toolkit并非仅做推理的强制安装项，驱动需满足相应PyTorch CUDA构建要求；重新编译第三方原生扩展需对应编译器。

冻结基座准备使用同一ESM环境，先安装上文CUDA12.4源torch2.6.0，再安装esm2.txt。转换本身只需CPU，严格复现原fair-esm/Transformers/safetensors版本；首次模型准备命令为`python tools/prepare_base.py`。


所有命令从解压后的发布根目录执行，输出目录需为空。不同代码块是独立示例，重复运行应换用新的输出目录。

若系统 Python 3.10 缺少 ensurepip（本轮服务器实际情况），不需要借用旧虚拟环境，也不需要修改系统；先在发布目录建立独立工作区，使用官方 pip 引导程序：

```bash
mkdir -p .runtime
curl --fail --location https://bootstrap.pypa.io/get-pip.py -o .runtime/get-pip.py
python3.10 -m venv --without-pip .env-post
.env-post/bin/python .runtime/get-pip.py pip==25.2
python3.10 -m venv --without-pip .env-esm2
.env-esm2/bin/python .runtime/get-pip.py pip==25.2
```

然后激活相应环境，从其上方安装代码块中的 pip install 步骤继续。Python必须为3.10；本轮使用系统3.10.12。下载脚本的本轮URL和SHA记录于验收日志 pip_bootstrap.json。该引导不改变科研依赖版本。
本轮服务器实际独立环境路径：`/mnt/ws/teeth/release_tests/IL24_20261001_183025/env-post` 与 `/mnt/ws/teeth/release_tests/IL24_20261001_183025/env-esm2`，它们不随ZIP交付；用户按上方步骤创建自己的环境即可。
