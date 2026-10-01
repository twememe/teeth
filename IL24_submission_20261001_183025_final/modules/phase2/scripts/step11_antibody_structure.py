#!/usr/bin/env python3
"""
Step 11: 抗体结构预测 - ABodyBuilder2
======================================
预测 IA6-13-8 抗体的 VH/VL Fv 三维结构

输入: VH_VL.fasta
输出: IA6-13-8_Fv.pdb, abodybuilder2_error_annotations.csv
"""

import os
import sys
import traceback
from pathlib import Path

try:
    import ImmuneBuilder
    HAS_IMMUNEBUILDER = True
except ImportError:
    HAS_IMMUNEBUILDER = False
    print("⚠️ ImmuneBuilder 未安装，将生成预测脚本供服务器使用")

import pandas as pd


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "results" / "11_antibody_structure"
INPUT_FASTA = DATA_DIR / "VH_VL_clean.fasta"

# 确保输出目录存在
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 序列解析
# ============================================================

def parse_fasta(fasta_path: Path) -> dict:
    """解析 FASTA 文件，提取 VH 和 VL 序列"""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def identify_chain(sequences: dict) -> tuple:
    """识别 VH 和 VL 序列"""
    vh_seq = None
    vl_seq = None
    vh_name = None
    vl_name = None
    
    for name, seq in sequences.items():
        # VH 通常以 V 开头的框架区开始
        if seq.startswith('EVQL') or seq.startswith('QVQL') or seq.startswith('IVQL'):
            vh_seq = seq
            vh_name = name
        # VL (kappa) 通常以 DIVMT 或 DVVMT 开始
        elif seq.startswith('DIVMT') or seq.startswith('DVVMT'):
            vl_seq = seq
            vl_name = name
    
    if vh_seq is None:
        raise ValueError("未找到 VH 序列！")
    if vl_seq is None:
        raise ValueError("未找到 VL 序列！")
    
    return vh_name, vh_seq, vl_name, vl_seq


# ============================================================
# ABodyBuilder2 预测
# ============================================================

def predict_with_immune_builder(vh_seq: str, vl_seq: str, output_dir: Path):
    """使用 ABodyBuilder2 预测抗体结构"""
    print("=" * 60)
    print("Step 11: ABodyBuilder2 抗体结构预测")
    print("=" * 60)
    
    # 初始化模型
    print("📦 加载 ABodyBuilder2 模型...")
    model = ImmuneBuilder.ABodyBuilder2()
    
    # 预测
    print(f"🔮 预测 VH ({len(vh_seq)} aa) + VL ({len(vl_seq)} aa)...")
    # ImmuneBuilder expects a chain-keyed sequence dictionary.
    structure = model.predict({'H': vh_seq, 'L': vl_seq})
    
    # 保存结构
    output_dir = output_dir / "ensemble"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pdb = output_dir / "IA6-13-8_Fv_abodybuilder2.pdb"
    # Current ImmuneBuilder returns Antibody.save_all(), which writes PDB.
    structure.save_all(str(output_dir))
    generated = sorted(output_dir.glob('*.pdb'))
    if generated:
        generated[0].replace(output_pdb)
    print(f"✅ 结构已保存: {output_pdb}")
    
    # 提取置信度分数
    # ImmuneBuilder 的误差估计并非 pLDDT；未经细化且 B 列全零时不作置信度解释
    # 我们从 PDB 文件的 b-factor 列读取
    confidence_data = extract_confidence_from_pdb(output_pdb, vh_seq, vl_seq)
    
    return output_pdb, confidence_data


def extract_confidence_from_pdb(pdb_path: Path, vh_seq: str, vl_seq: str) -> pd.DataFrame:
    """从 PDB 文件提取置信度分数"""
    data = []
    current_chain = None
    residue_num = 0
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                b_factor = float(line[60:66].strip())
                
                # 只记录 CA 原子
                if line[12:16].strip() == 'CA':
                    data.append({
                        'chain': chain,
                        'residue': res_num,
                        'residue_name': res_name,
                        'pdb_b_column': b_factor,
                        'metric': 'ABodyBuilder2_predicted_error_if_written',
                        'interpretation': 'not_pLDDT; zero_in_unrefined_files_is_not_confidence'
                    })
    
    df = pd.DataFrame(data)
    
    # 保存置信度表
    confidence_file = pdb_path.parent / "abodybuilder2_error_annotations.csv"
    df.to_csv(confidence_file, index=False)
    print(f"📊 置信度分数已保存: {confidence_file}")
    
    return df


# ============================================================
# 备选：生成服务器脚本
# ============================================================

def generate_server_script(vh_seq: str, vl_seq: str, output_dir: Path):
    """生成可在服务器上运行的脚本"""
    
    script_content = f'''#!/bin/bash
# ABodyBuilder2 抗体结构预测脚本
# ============================================================
# 在服务器上运行此脚本前，请先安装依赖：
#   conda create -n immune python=3.10
#   conda activate immune
#   conda install -c conda-forge pytorch openmm pdbfixer
#   conda install -c bioconda anarci
#   pip install ImmuneBuilder

# 创建输出目录
mkdir -p {output_dir}

# Python 预测脚本
python << 'EOF'
import ImmuneBuilder

# VH 序列 ({len(vh_seq)} aa)
VH_SEQ = "{vh_seq}"

# VL 序列 ({len(vl_seq)} aa)  
VL_SEQ = "{vl_seq}"

print("Loading ABodyBuilder2 model...")
model = ImmuneBuilder.ABodyBuilder2()

print("Predicting structure...")
structure = model.predict(heavy=VH_SEQ, light=VL_SEQ)

output_path = "{output_dir}/IA6-13-8_Fv.pdb"
structure.save_pdb(output_path)
print(f"Structure saved: {{output_path}}")
EOF
'''
    
    script_file = output_dir / "run_abodybuilder2.sh"
    with open(script_file, 'w') as f:
        f.write(script_content)
    os.chmod(script_file, 0o755)
    print(f"📝 服务器脚本已生成: {script_file}")


def generate_docker_script(vh_seq: str, vl_seq: str, output_dir: Path):
    """生成 Docker 运行脚本"""
    
    script_content = f'''#!/bin/bash
# ABodyBuilder2 Docker 运行脚本
# ============================================================
# 确保 Docker 已安装，NVIDIA Container Toolkit 已配置

# 创建数据目录
mkdir -p /tmp/abodybuilder2_data

# 写入 FASTA 文件
cat > /tmp/abodybuilder2_data/VH_VL.fasta << 'FASTA'
>IA6-13-8_VH
{vh_seq}
>IA6-13-8_VL
{vl_seq}
FASTA

# 运行 Docker 容器
docker run -it \\
  -v /tmp/abodybuilder2_data:/data \\
  oxpig/abodybuilder2:latest \\
  ABodyBuilder2 \\
  --fasta_file /data/VH_VL.fasta \\
  --output_dir /data/output

# 复制结果到当前目录
cp -r /tmp/abodybuilder2_data/output/* {output_dir}/

echo "✅ 结果已保存到: {output_dir}"
'''
    
    script_file = output_dir / "run_docker.sh"
    with open(script_file, 'w') as f:
        f.write(script_content)
    os.chmod(script_file, 0o755)
    print(f"🐳 Docker 脚本已生成: {script_file}")


# ============================================================
# 主函数
# ============================================================

def main():
    print("=" * 60)
    print("Step 11: 抗体结构预测")
    print("=" * 60)
    
    # 检查输入文件
    if not INPUT_FASTA.exists():
        print(f"❌ 错误: 找不到输入文件 {INPUT_FASTA}")
        print("请先运行 Step 10 或确保 VH_VL.fasta 存在")
        sys.exit(1)
    
    # 解析序列
    print(f"📂 读取序列: {INPUT_FASTA}")
    sequences = parse_fasta(INPUT_FASTA)
    print(f"   找到 {len(sequences)} 个序列:")
    for name, seq in sequences.items():
        print(f"   - {name}: {len(seq)} aa")
    
    # 识别 VH/VL
    vh_name, vh_seq, vl_name, vl_seq = identify_chain(sequences)
    print(f"\n🔍 VH: {vh_name} ({len(vh_seq)} aa)")
    print(f"🔍 VL: {vl_name} ({len(vl_seq)} aa)")
    
    # 保存解析后的序列
    seq_info_file = OUTPUT_DIR / "sequence_info.txt"
    with open(seq_info_file, 'w') as f:
        f.write(f"VH: {vh_name}\\n")
        f.write(f"Sequence ({len(vh_seq)} aa):\\n{vh_seq}\\n\\n")
        f.write(f"VL: {vl_name}\\n")
        f.write(f"Sequence ({len(vl_seq)} aa):\\n{vl_seq}\\n")
    print(f"📄 序列信息已保存: {seq_info_file}")
    
    if HAS_IMMUNEBUILDER:
        # 直接运行预测
        try:
            output_pdb, confidence_df = predict_with_immune_builder(
                vh_seq, vl_seq, OUTPUT_DIR
            )
            
            # 输出统计
            print("\\n📊 置信度统计:")
            print(confidence_df[['pdb_b_column','interpretation']].describe(include='all'))
            
        except Exception as e:
            print(f"⚠️ 预测失败: {e!r}")
            traceback.print_exc()
            print("将生成服务器脚本...")
            generate_server_script(vh_seq, vl_seq, OUTPUT_DIR)
            generate_docker_script(vh_seq, vl_seq, OUTPUT_DIR)
    else:
        print("\\n⚠️ ImmuneBuilder 未安装")
        print("生成服务器脚本...")
        generate_server_script(vh_seq, vl_seq, OUTPUT_DIR)
        generate_docker_script(vh_seq, vl_seq, OUTPUT_DIR)
    
    print("\\n" + "=" * 60)
    print("Step 11 完成!")
    print("=" * 60)
    print(f"输出目录: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
