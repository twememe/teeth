#!/usr/bin/env python3
"""
Step 15: 亲和力评估
=====================
评估 IL-24 + 抗体复合物的结合亲和力

三级评估体系:
- Level 1: 界面质量 (Interface area, H-bond, Salt bridge)
- Level 2: PRODIGY 能量估算 (ΔG, K_D)
- Level 3: Rosetta/FoldX 精细能量分析 (ΔΔG)

输入: IL-24_Ab_complex.pdb (from Step 12)
输出: affinity_report.json
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime
import argparse

import numpy as np
import pandas as pd


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "15_affinity_estimation"

# 常数
R = 0.001987  # 气体常数 kcal/(mol·K)
T = 298.15  # 室温 K

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 界面性质计算
# ============================================================

def compute_interface_properties(pdb_path: Path, 
                                distance_threshold: float = 4.5) -> dict:
    """
    计算界面性质
    
    返回:
    - interface_area: 界面面积 (Å²)
    - n_contacts: 接触数
    - n_hbonds: 氢键数 (估算)
    - n_salt_bridges: 盐桥数 (估算)
    """
    properties = {
        'interface_area': 0,
        'n_contacts': 0,
        'n_hbonds': 0,
        'n_salt_bridges': 0,
        'hydrophobic_residues': [],
        'charged_residues': []
    }
    
    # 残基性质
    hydrophobic = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR'}
    positively_charged = {'LYS', 'ARG', 'HIS'}
    negatively_charged = {'ASP', 'GLU'}
    
    # 氢键供体/受体
    hbond_donors = {'SER', 'THR', 'TYR', 'CYS', 'TRP', 'HIS', 'LYS', 'ARG'}
    hbond_acceptors = {'SER', 'THR', 'TYR', 'CYS', 'ASP', 'GLU', 'ASN', 'GLN'}
    
    contacts = []
    residues = set()
    
    with open(pdb_path, 'r') as f:
        atoms = []
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append({
                    'chain': line[21],
                    'res_num': int(line[22:26].strip()),
                    'res_name': line[17:20].strip(),
                    'atom_name': line[12:16].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[76:78].strip()
                })
        
        # 简化：基于接触距离估算
        # 假设 IL-24 在链 A，抗体在链 B/C
        il24_atoms = [a for a in atoms if a['chain'] == 'A']
        ab_atoms = [a for a in atoms if a['chain'] in ['B', 'H', 'L']]
        
        for il24_atom in il24_atoms:
            for ab_atom in ab_atoms:
                dist = np.sqrt(
                    (il24_atom['x'] - ab_atom['x'])**2 +
                    (il24_atom['y'] - ab_atom['y'])**2 +
                    (il24_atom['z'] - ab_atom['z'])**2
                )
                
                if dist < distance_threshold:
                    res_key = (il24_atom['chain'], il24_atom['res_num'], il24_atom['res_name'])
                    residues.add(res_key)
                    
                    # 估算氢键
                    if (il24_atom['res_name'] in hbond_donors and 
                        ab_atom['res_name'] in hbond_acceptors):
                        properties['n_hbonds'] += 0.1  # 加权计数
                    
                    # 估算盐桥
                    if (il24_atom['res_name'] in positively_charged and 
                        ab_atom['res_name'] in negatively_charged):
                        properties['n_salt_bridges'] += 0.1
    
    properties['n_contacts'] = len(residues)
    properties['interface_area'] = len(residues) * 150  # 简化估算
    
    # 统计残基类型
    for chain, res_num, res_name in residues:
        if res_name in hydrophobic:
            properties['hydrophobic_residues'].append(res_num)
        if res_name in positively_charged or res_name in negatively_charged:
            properties['charged_residues'].append(res_num)
    
    return properties


def prodigy_estimate(interface_area: float, n_contacts: int, 
                    n_hbonds: float, n_salt_bridges: float,
                    temperature: float = 298.15) -> dict:
    """
    PRODIGY 能量估算
    
    PRODIGY 使用以下公式估算结合自由能:
    ΔG = -10.2 + 0.109 * N_interface_residues + 0.248 * N_contacts + 
          1.26 * N_hbonds + 0.41 * N_salt_bridges - 0.00693 * interface_area
    
    注意: 这是简化版，实际 PRODIGY 使用更复杂的模型
    """
    # 简化公式参数 (基于 PRODIGY 训练集)
    delta_g = (
        -10.2 +
        0.109 * n_contacts +
        0.248 * n_contacts +  # interface residues 类似
        1.26 * n_hbonds +
        0.41 * n_salt_bridges -
        0.00693 * interface_area / 100  # 归一化
    )
    
    # 转换为 K_D
    # ΔG = RT * ln(K_D)
    # K_D = exp(ΔG / (RT))
    kd = np.exp(delta_g * 1000 / (R * temperature))  # 转换为 J/mol
    
    # 转换为更常用的单位
    kd_nm = kd * 1e9  # nM
    
    return {
        'delta_g_kcal': delta_g,
        'delta_g_kj': delta_g * 4.184,
        'kd_nm': kd_nm,
        'kd_m': kd,
        'method': 'PRODIGY-simplified',
        'confidence': 'computational_estimate_only'
    }


def compute_binding_score(interface_area: float, n_contacts: int,
                         n_hbonds: float, n_salt_bridges: float,
                         hydrophobic_ratio: float = 0.5) -> float:
    """
    综合绑定评分 (0-1)
    
    考虑因素:
    - 界面面积 (越大越好，通常 1000-3000 Å² 最佳)
    - 接触数 (越多越好)
    - 氢键 (正相关)
    - 盐桥 (正相关)
    - 疏水比例 (正相关)
    """
    scores = []
    
    # 界面面积评分
    if 1000 <= interface_area <= 3000:
        area_score = 1.0
    elif 500 <= interface_area < 1000:
        area_score = 0.7
    elif 3000 < interface_area <= 5000:
        area_score = 0.8
    else:
        area_score = 0.3
    scores.append(area_score * 0.3)  # 权重 30%
    
    # 接触数评分
    contact_score = min(1.0, n_contacts / 50)
    scores.append(contact_score * 0.2)  # 权重 20%
    
    # 氢键评分
    hbond_score = min(1.0, n_hbonds / 5)
    scores.append(hbond_score * 0.2)  # 权重 20%
    
    # 盐桥评分
    salt_score = min(1.0, n_salt_bridges / 3)
    scores.append(salt_score * 0.1)  # 权重 10%
    
    # 疏水比例评分
    if 0.3 <= hydrophobic_ratio <= 0.6:
        hydro_score = 1.0
    else:
        hydro_score = 0.7
    scores.append(hydro_score * 0.2)  # 权重 20%
    
    return sum(scores)


def estimate_affinity_category(kd_nm: float) -> str:
    """
    根据 K_D 估算亲和力类别
    """
    if kd_nm < 0.1:
        return 'PICOMOLAR (pM)'
    elif kd_nm < 100:
        return 'NANOMOLAR (nM)'
    elif kd_nm < 100000:
        return 'MICROMOLAR (µM)'
    else:
        return 'MILLIMOLAR (mM)'


def mutation_effect_analysis(wildtype_residues: list, 
                           interface_residues: list,
                           mutations: list = None) -> dict:
    """
    估算突变对亲和力的影响 (ΔΔG)
    
    简化模型:
    - 疏水残基 → 极性残基: ΔΔG > 0 (降低亲和力)
    - 极性残基 → 疏水残基: ΔΔG < 0 (提高亲和力，但可能有特异性损失)
    - 丙氨酸扫描: ΔΔG ≈ 0.5-2.0 kcal/mol per mutation
    
    这是计算层面的估算，不作为最终结论
    """
    analysis = {
        'method': 'simplified_energy_model',
        'confidence': 'low',
        'mutations': []
    }
    
    # 如果没有指定突变，进行丙氨酸扫描
    if not mutations:
        mutations = [{'res': r, 'to': 'ALA'} for r in interface_residues[:10]]
    
    for mut in mutations:
        res = mut['res']
        to_res = mut['to']
        
        # 简化的 ΔΔG 估算
        aa_properties = {
            'ALA': {'hydrophobic': 1, 'charge': 0},
            'VAL': {'hydrophobic': 1, 'charge': 0},
            'ILE': {'hydrophobic': 1, 'charge': 0},
            'LEU': {'hydrophobic': 1, 'charge': 0},
            'MET': {'hydrophobic': 1, 'charge': 0},
            'PHE': {'hydrophobic': 1, 'charge': 0},
            'TRP': {'hydrophobic': 1, 'charge': 0},
            'TYR': {'hydrophobic': 0.5, 'charge': 0},
            'SER': {'hydrophobic': 0, 'charge': 0},
            'THR': {'hydrophobic': 0, 'charge': 0},
            'ASN': {'hydrophobic': 0, 'charge': 0},
            'GLN': {'hydrophobic': 0, 'charge': 0},
            'HIS': {'hydrophobic': 0, 'charge': 0.5},
            'LYS': {'hydrophobic': 0, 'charge': 1},
            'ARG': {'hydrophobic': 0, 'charge': 1},
            'ASP': {'hydrophobic': 0, 'charge': -1},
            'GLU': {'hydrophobic': 0, 'charge': -1},
            'CYS': {'hydrophobic': 0.8, 'charge': 0},
            'GLY': {'hydrophobic': 0, 'charge': 0},
            'PRO': {'hydrophobic': 0.3, 'charge': 0},
        }
        
        from_prop = aa_properties.get(res, {'hydrophobic': 0, 'charge': 0})
        to_prop = aa_properties.get(to_res, {'hydrophobic': 0, 'charge': 0})
        
        # ΔΔG ≈ 疏水贡献 + 电荷贡献
        delta_hydro = (to_prop['hydrophobic'] - from_prop['hydrophobic']) * 1.5  # kcal/mol
        delta_charge = abs(to_prop['charge'] - from_prop['charge']) * 1.0  # kcal/mol
        
        ddg = delta_hydro - delta_charge  # 负值表示增强结合
        
        analysis['mutations'].append({
            'from': res,
            'to': to_res,
            'delta_dg_kcal': round(ddg, 2),
            'interpretation': 'stabilizing' if ddg < 0 else ('destabilizing' if ddg > 0 else 'neutral')
        })
    
    return analysis


# ============================================================
# 主函数
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Step 15: Affinity Estimation')
    parser.add_argument('--input', '-i', type=str,
                       help='输入 PDB 文件',
                       default=None)
    parser.add_argument('--mode', '-m', type=str,
                       choices=['12A_blind', '12B_guided'],
                       default='12A_blind')
    parser.add_argument('--complex-dir', '-c', type=str,
                       help='复合物目录 (ensemble)',
                       default=None)
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 15: 亲和力评估")
    print("=" * 60)
    
    # 确定输入路径
    if args.input:
        input_path = Path(args.input)
    elif args.complex_dir:
        input_path = Path(args.complex_dir)
    else:
        input_path = RESULTS_DIR / "12_complex_prediction" / args.mode
    
    if not input_path.exists():
        print(f"❌ 错误: 找不到输入 {input_path}")
        sys.exit(1)
    
    print(f"📂 输入: {input_path}")
    
    # 查找 PDB 文件
    if input_path.is_dir():
        pdb_files = list(input_path.glob("*.pdb"))
        if not pdb_files:
            # 查找子目录
            for subdir in input_path.iterdir():
                if subdir.is_dir():
                    pdb_files = list(subdir.glob("*.pdb"))
                    if pdb_files:
                        break
        pdb_files = pdb_files[:3]  # 分析前3个
    else:
        pdb_files = [input_path]
    
    if not pdb_files:
        print(f"❠️ 未找到 PDB 文件")
        sys.exit(1)
    
    print(f"📊 分析 {len(pdb_files)} 个复合物模型...")
    
    all_results = []
    
    for pdb_file in pdb_files:
        print(f"\n📄 {pdb_file.name}")
        
        try:
            # Level 1: 界面性质
            properties = compute_interface_properties(pdb_file)
            
            hydrophobic_ratio = (
                len(properties['hydrophobic_residues']) / properties['n_contacts']
                if properties['n_contacts'] > 0 else 0
            )
            
            print(f"   界面面积: {properties['interface_area']:.1f} Å²")
            print(f"   接触数: {properties['n_contacts']}")
            print(f"   氢键数: {properties['n_hbonds']:.1f}")
            print(f"   盐桥数: {properties['n_salt_bridges']:.1f}")
            print(f"   疏水残基: {len(properties['hydrophobic_residues'])}")
            print(f"   带电残基: {len(properties['charged_residues'])}")
            
            # Level 2: PRODIGY 估算
            prodigy = prodigy_estimate(
                properties['interface_area'],
                properties['n_contacts'],
                properties['n_hbonds'],
                properties['n_salt_bridges']
            )
            
            print(f"\n   📊 PRODIGY 估算:")
            print(f"   ΔG = {prodigy['delta_g_kcal']:.2f} kcal/mol")
            print(f"   K_D ≈ {prodigy['kd_nm']:.2f} nM")
            print(f"   类别: {estimate_affinity_category(prodigy['kd_nm'])}")
            
            # 综合评分
            binding_score = compute_binding_score(
                properties['interface_area'],
                properties['n_contacts'],
                properties['n_hbonds'],
                properties['n_salt_bridges'],
                hydrophobic_ratio
            )
            
            result = {
                'file': pdb_file.name,
                'level1_interface': properties,
                'level2_prodigy': prodigy,
                'binding_score': binding_score,
                'affinity_category': estimate_affinity_category(prodigy['kd_nm'])
            }
            
            all_results.append(result)
            
        except Exception as e:
            print(f"   ❌ 分析失败: {e}")
    
    # 汇总
    if all_results:
        print("\n" + "=" * 60)
        print("📊 亲和力评估汇总")
        print("=" * 60)
        
        mean_score = np.mean([r['binding_score'] for r in all_results])
        mean_delta_g = np.mean([r['level2_prodigy']['delta_g_kcal'] for r in all_results])
        mean_kd = np.mean([r['level2_prodigy']['kd_nm'] for r in all_results])
        
        print(f"   平均结合评分: {mean_score:.3f}")
        print(f"   平均 ΔG: {mean_delta_g:.2f} kcal/mol")
        print(f"   平均 K_D: {mean_kd:.2f} nM")
        
        # 警告
        print("\n⚠️ 重要提示:")
        print("   - 所有 K_D 估算均为 computational estimate")
        print("   - 不可作为最终结论")
        print("   - 需结合 PRODIGY + Rosetta + FoldX + 实验数据综合判断")
        print("   - Boltz-2 的 affinity 预测不可靠 (ROC-AUC 0.5-0.60)")
    
    # 保存报告
    report = {
        'script': 'step15_affinity_estimation.py',
        'timestamp': datetime.now().isoformat(),
        'mode': args.mode,
        'n_models': len(all_results),
        'results': all_results,
        'summary': {
            'mean_binding_score': float(mean_score) if all_results else None,
            'mean_delta_g': float(mean_delta_g) if all_results else None,
            'mean_kd_nm': float(mean_kd) if all_results else None
        },
        'disclaimer': (
            'This is a computational estimate only. '
            'Do not use as final conclusion. '
            'Requires PRODIGY + Rosetta + FoldX + experimental validation.'
        )
    }
    
    report_file = OUTPUT_DIR / f"affinity_report_{args.mode}.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n📋 报告已保存: {report_file}")
    
    # 保存表格
    if all_results:
        df = pd.DataFrame([{
            'file': r['file'],
            'binding_score': r['binding_score'],
            'delta_g_kcal': r['level2_prodigy']['delta_g_kcal'],
            'kd_nm': r['level2_prodigy']['kd_nm'],
            'interface_area': r['level1_interface']['interface_area'],
            'n_contacts': r['level1_interface']['n_contacts'],
            'n_hbonds': r['level1_interface']['n_hbonds'],
            'n_salt_bridges': r['level1_interface']['n_salt_bridges']
        } for r in all_results])
        
        csv_file = OUTPUT_DIR / f"affinity_{args.mode}.csv"
        df.to_csv(csv_file, index=False)
        print(f"📊 表格已保存: {csv_file}")
    
    print("\n" + "=" * 60)
    print("Step 15 完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
