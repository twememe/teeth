#!/usr/bin/env python3
"""
Step 12: IL-24 + 抗体复合物预测
=================================
包含两种模式:
- 12A: 纯盲复合物预测 (无先验)
- 12B: Phase1-Guided 预测 (使用 Phase 1 热点作为软约束)

输入: IL-24 结构 (PDB), VH/VL Fv 结构 (PDB)
输出: IL-24_Ab_complex.pdb (ensemble)
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime

try:
    from openfold.np import protein
    HAS_OPENFOLD = True
except ImportError:
    HAS_OPENFOLD = False


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "results" / "12_complex_prediction"
IL24_PDB = DATA_DIR / "AF-Q925S4-F1-model_v6.pdb"  # IL-24 AlphaFold 结构
ANTIBODY_FV = OUTPUT_DIR / "11_antibody_structure" / "IA6-13-8_Fv.pdb"

# Phase 1 热点区域 (用于 Guided 预测)
PHASE1_HOTSPOTS = [
    {"region": "C1-B", "start": 52, "end": 66, "confidence": "VERY_HIGH"},
    {"region": "C2-B", "start": 126, "end": 140, "confidence": "VERY_HIGH"},
    {"region": "C3-B", "start": 60, "end": 74, "confidence": "VERY_HIGH"},
    {"region": "C4-B", "start": 143, "end": 157, "confidence": "VERY_HIGH"},
    {"region": "C5-B", "start": 117, "end": 131, "confidence": "VERY_HIGH"},
]

# 软约束权重
SOFT_RESTRAINT_WEIGHT = 0.3  # 建议值 0.2-0.5

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 辅助函数
# ============================================================

def read_fasta(fasta_path: Path) -> dict:
    """读取 FASTA 文件"""
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


def parse_pdb(pdb_path: Path) -> dict:
    """解析 PDB 文件，提取链和残基信息"""
    chains = {}
    current_chain = None
    current_residues = []
    current_res_num = None
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                
                if chain != current_chain:
                    if current_chain is not None:
                        chains[current_chain] = current_residues
                    current_chain = chain
                    current_residues = []
                
                current_residues.append({
                    'res_num': res_num,
                    'res_name': res_name,
                    'x': x, 'y': y, 'z': z
                })
    
    if current_chain is not None:
        chains[current_chain] = current_residues
    
    return chains


def generate_phase1_constraints(hotspots: list, il24_chain: str = 'A') -> dict:
    """
    将 Phase 1 热点转换为模型约束
    
    返回格式 (可用于 Chai-1 或 AF3):
    - residue_id: 约束类型和权重
    """
    constraints = {}
    
    for hotspot in hotspots:
        region = hotspot['region']
        start = hotspot['start']
        end = hotspot['end']
        confidence = hotspot['confidence']
        
        # 权重基于置信度
        weight_map = {
            'VERY_HIGH': 1.0,
            'HIGH': 0.8,
            'MEDIUM': 0.5,
            'LOW': 0.3
        }
        weight = weight_map.get(confidence, 0.5) * SOFT_RESTRAINT_WEIGHT
        
        for res in range(start, end + 1):
            constraints[f"{il24_chain}:{res}"] = {
                'type': 'surface_exposure',
                'weight': weight,
                'region': region,
                'source': 'Phase1'
            }
    
    return constraints


def generate_af3_input(il24_pdb: Path, antibody_fv: Path, output_dir: Path, 
                       use_constraints: bool = False, constraints: dict = None):
    """
    生成 AlphaFold3 输入文件
    
    AF3 使用 JSON 格式定义复合物:
    {
        "name": "IL24_Ab_complex",
        "sequences": [
            {"protein": {"name": "IL24", "sequence": "XXX"}},
            {"protein": {"name": "VH", "sequence": "XXX"}},
            {"protein": {"name": "VL", "sequence": "XXX"}}
        ],
        " restraints": [...]  # 可选
    }
    """
    # 读取 IL-24 序列
    with open(il24_pdb, 'r') as f:
        il24_seq = []
        current_res = None
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                res = int(line[22:26].strip())
                if res != current_res:
                    il24_seq.append(line[17:20].strip())
                    current_res = res
        il24_sequence = ''.join([AA3_TO_1.get(aa, 'X') for aa in il24_seq])
    
    # 从 antibody PDB 提取序列
    if antibody_fv.exists():
        ab_chains = parse_pdb(antibody_fv)
        vh_seq = ''.join([AA3_TO_1.get(r['res_name'], 'X') 
                          for r in ab_chains.get('H', [])])
        vl_seq = ''.join([AA3_TO_1.get(r['res_name'], 'X') 
                          for r in ab_chains.get('L', [])])
    else:
        # 使用原始 FASTA
        fasta = read_fasta(DATA_DIR / "VH_VL.fasta")
        vh_seq = [v for k, v in fasta.items() if 'VH' in k][0]
        vl_seq = [v for k, v in fasta.items() if 'VL' in k][0]
    
    # 构建 AF3 JSON
    af3_input = {
        "name": "IL24_IA6-13-8_complex",
        "modelSeeds": list(range(1, 21)),  # 20 seeds for ensemble
        "sequences": [
            {
                "protein": {
                    "name": "IL24",
                    "sequence": il24_sequence
                }
            },
            {
                "protein": {
                    "name": "VH",
                    "sequence": vh_seq
                }
            },
            {
                "protein": {
                    "name": "VL",
                    "sequence": vl_seq
                }
            }
        ]
    }
    
    # 添加约束（如果启用）
    if use_constraints and constraints:
        af3_input["dialect"] = "alphafold"
        af3_input["version"] = 1
        af3_input["restraints"] = []
        
        for key, value in constraints.items():
            chain, res_num = key.split(':')
            af3_input["restraints"].append({
                "type": "custom_polymer_restraint",
                "weight": value['weight'],
                "chains": [chain],
                "residueIndex": [int(res_num)],
                "feature": "surface_exposure",
                "flatness": 0.5,
                "label": f"Phase1_{value['region']}"
            })
    
    return af3_input


def generate_chai1_constraints(hotspots: list, il24_chain: str = 'A') -> str:
    """
    生成 Chai-1 的 MSA 约束格式
    
    Chai-1 支持自定义 MSA，可以通过约束增强 Phase 1 热点区域
    """
    constraint_lines = []
    
    for hotspot in hotspots:
        region = hotspot['region']
        start = hotspot['start']
        end = hotspot['end']
        weight = SOFT_RESTRAINT_WEIGHT * (
            1.0 if hotspot['confidence'] == 'VERY_HIGH' else 
            0.8 if hotspot['confidence'] == 'HIGH' else 0.5
        )
        
        constraint_lines.append(
            f"# {region}: {start}-{end} (weight={weight:.2f})"
        )
        constraint_lines.append(
            f"A:{start}-{end} SASA > 0.25"
        )
    
    return '\n'.join(constraint_lines)


# AA 3-letter to 1-letter conversion
AA3_TO_1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


# ============================================================
# 生成服务器执行脚本
# ============================================================

def generate_af3_server_script(af3_input: dict, mode: str, output_dir: Path):
    """生成 AF3 服务器执行脚本"""
    
    mode_label = "12A_blind" if mode == "blind" else "12B_guided"
    script_path = output_dir / f"run_af3_{mode_label}.sh"
    
    # 将输入 JSON 保存
    json_path = output_dir / f"af3_input_{mode_label}.json"
    with open(json_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    script_content = f'''#!/bin/bash
# AlphaFold3 复合物预测脚本 (模式: {mode_label.upper()})
# ============================================================

# 配置
MODE="{mode_label}"
OUTPUT_DIR="{output_dir}"
JSON_INPUT="{json_path}"

# Step 1: 准备输入 (已在本地生成)
echo "📂 输入文件: ${{JSON_INPUT}}"

# Step 2: 上传到 AlphaFold3 服务器
# 方法1: 使用 AlphaFold3 Web API (如果可用)
# 方法2: 使用本地 AlphaFold3 安装
# 方法3: 使用 Google Colab (备用)

# 检查是否有 AF3 本地安装
if command -v af3_predict &> /dev/null; then
    echo "🔮 使用本地 AF3..."
    af3_predict \\
        --json_input "${{JSON_INPUT}}" \\
        --output_dir "${{OUTPUT_DIR}}/af3_output_{mode_label}" \\
        --num_seeds 20
    
elif command -v python &> /dev/null; then
    echo "🔮 使用 Python API (AlphaFold3 或 Boltz-1)..."
    
    # 示例: 使用 Boltz-1
    python << 'AF3_EOF'
import json
import sys

try:
    from boltz import BoltzModel
    
    with open("{json_path}") as f:
        input_data = json.load(f)
    
    print("Loading Boltz-1 model...")
    model = BoltzModel()
    
    print("Running prediction (20 seeds)...")
    results = model.predict(
        sequences=input_data["sequences"],
        restraints=input_data.get("restraints", []),
        num_seeds=20
    )
    
    # 保存每个 seed 的结果
    output_dir = "{output_dir}/af3_output_{mode_label}"
    os.makedirs(output_dir, exist_ok=True)
    
    for i, result in enumerate(results):
        result.save_pdb(f"{{output_dir}}/rank_{{i}}.pdb")
    
    print(f"✅ Done! Results saved to: {{output_dir}}")
    
except ImportError as e:
    print(f"⚠️ 需要安装: pip install boltz")
    print(f"⚠️ 错误: {{e}}")
AF3_EOF

else
    echo "⚠️ 未找到 AF3 或 Boltz-1"
    echo "请手动上传 {json_path} 到 AlphaFold3 服务器"
    echo "服务器地址: https://alphafold.ebi.ac.uk/submit"
fi

# Step 3: 下载结果后，运行本地的后处理脚本
echo ""
echo "📋 下一步: 运行 step13_contact_mapping.py"
echo "   python scripts/step13_contact_mapping.py --mode {mode_label}"
'''
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"📝 AF3 执行脚本: {script_path}")


def generate_colab_notebook(af3_input: dict, mode: str, output_dir: Path):
    """生成 Google Colab notebook"""
    
    mode_label = "12A_blind" if mode == "blind" else "12B_guided"
    nb_path = output_dir / f"af3_{mode_label}.ipynb"
    
    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.10.0"
            }
        },
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [f"# AlphaFold3 复合物预测 - {mode_label.upper()}\n\n"]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# 安装依赖\n",
                    "!pip install -q boltz alphafold3\n",
                    "import json\n",
                    "\n",
                    f"# 输入数据 (mode: {mode_label})\n",
                    f"input_data = {json.dumps(af3_input, indent=2)}\n",
                    "\n",
                    "print(f'Sequences: {{len(input_data[\"sequences\"])}} chains')\n",
                    "if 'restraints' in input_data:\n",
                    "    print(f'Restraints: {{len(input_data[\"restraints\"])}}')"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# 运行 Boltz-1 预测\n",
                    "from boltz import BoltzModel\n",
                    "import os\n",
                    "\n",
                    "print('Loading Boltz-1 model...')\n",
                    "model = BoltzModel()\n",
                    "\n",
                    "print('Running prediction (20 seeds)...')\n",
                    "results = model.predict(\n",
                    "    sequences=input_data['sequences'],\n",
                    "    restraints=input_data.get('restraints', []),\n",
                    "    num_seeds=20\n",
                    ")\n",
                    "\n",
                    f"# 保存结果\n",
                    f"output_dir = '{output_dir}/af3_output_{mode_label}'\n",
                    "os.makedirs(output_dir, exist_ok=True)\n",
                    "for i, result in enumerate(results):\n",
                    "    result.save_pdb(f'{output_dir}/rank_{{i}}.pdb')\n",
                    "\n",
                    "print(f'✅ Done! {{len(results)}} models saved')"
                ]
            }
        ]
    }
    
    with open(nb_path, 'w') as f:
        json.dump(notebook, f, indent=2)
    print(f"📓 Colab notebook: {nb_path}")


def generate_provenance(mode: str, output_dir: Path):
    """生成 provenance.json"""
    provenance = {
        "script": "step12_complex_prediction.py",
        "mode": mode,
        "timestamp": datetime.now().isoformat(),
        "inputs": {
            "il24_structure": str(IL24_PDB),
            "antibody_fv": str(ANTIBODY_FV),
            "phase1_hotspots": PHASE1_HOTSPOTS
        },
        "parameters": {
            "num_seeds": 20,
            "soft_restraint_weight": SOFT_RESTRAINT_WEIGHT
        },
        "outputs": {
            "af3_input": f"af3_input_{mode}.json",
            "scripts": [
                f"run_af3_{mode}.sh",
                f"af3_{mode}.ipynb"
            ]
        }
    }
    
    prov_path = output_dir / f"provenance_{mode}.json"
    with open(prov_path, 'w') as f:
        json.dump(provenance, f, indent=2)
    print(f"📋 Provenance: {prov_path}")


# ============================================================
# 主函数
# ============================================================

def main():
    print("=" * 60)
    print("Step 12: IL-24 + 抗体复合物预测")
    print("=" * 60)
    
    # 检查 IL-24 结构
    if not IL24_PDB.exists():
        print(f"❌ 错误: 找不到 IL-24 结构 {IL24_PDB}")
        sys.exit(1)
    
    print(f"📂 IL-24 结构: {IL24_PDB}")
    print(f"📂 抗体 Fv 结构: {ANTIBODY_FV}")
    print()
    
    # Phase 1 热点摘要
    print("🎯 Phase 1 热点区域:")
    for hs in PHASE1_HOTSPOTS:
        print(f"   {hs['region']}: {hs['start']}-{hs['end']} ({hs['confidence']})")
    print()
    
    # 生成两种模式的输入
    for mode in ["blind", "guided"]:
        print(f"\n{'='*60}")
        print(f"模式: {mode.upper()}")
        print('='*60)
        
        mode_label = "12A_blind" if mode == "blind" else "12B_guided"
        mode_dir = OUTPUT_DIR / mode_label
        mode_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成约束
        if mode == "guided":
            constraints = generate_phase1_constraints(PHASE1_HOTSPOTS)
            print(f"🔒 生成 Phase1-Guided 约束: {len(constraints)} 残基")
        else:
            constraints = None
            print("🔓 盲预测模式 (无先验信息)")
        
        # 生成 AF3 输入
        af3_input = generate_af3_input(
            IL24_PDB, ANTIBODY_FV, mode_dir,
            use_constraints=(mode == "guided"),
            constraints=constraints
        )
        
        # 生成执行脚本
        generate_af3_server_script(af3_input, mode, mode_dir)
        generate_colab_notebook(af3_input, mode, mode_dir)
        generate_provenance(mode, mode_dir)
        
        # 保存输入 JSON
        json_path = mode_dir / f"af3_input_{mode}.json"
        with open(json_path, 'w') as f:
            json.dump(af3_input, f, indent=2)
        
        print(f"\n📊 序列信息:")
        for seq in af3_input['sequences']:
            name = seq['protein']['name']
            seq_str = seq['protein']['sequence'][:50]
            length = len(seq['protein']['sequence'])
            print(f"   {name}: {length} aa (前50: {seq_str}...)")
    
    print("\n" + "=" * 60)
    print("Step 12 完成!")
    print("=" * 60)
    print(f"\n📁 输出目录: {OUTPUT_DIR}")
    print("\n下一步:")
    print("  1. 上传 phase2 目录到云服务器")
    print("  2. 在服务器上运行 run_af3_12A_blind.sh 或 run_af3_12B_guided.sh")
    print("  3. 下载结果后运行 step13_contact_mapping.py")


if __name__ == "__main__":
    main()
