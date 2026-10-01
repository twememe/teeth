#!/usr/bin/env python3
"""
Step 14: 复合物质量评估
=========================
评估 IL-24 + 抗体复合物的质量

评估指标:
- ipTM (interface Predicted TM-score)
- PAE (Predicted Aligned Error)
- Interface contacts
- CDR involvement
- Buried surface area
- Clashes
- Epitope consistency

输入: IL-24_Ab_complex.pdb (from Step 12)
输出: quality_report.json, metrics.csv
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import argparse

import numpy as np
import pandas as pd

try:
    from Bio import PDB
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("⚠️ 建议安装: pip install biopython")


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "14_quality_assessment"

# 质量阈值
CLASH_THRESHOLD = 0.4  # Å (两个原子的最小距离)
HIGH_CONFIDENCE_IPTM = 0.8
MEDIUM_CONFIDENCE_IPTM = 0.5

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# PDB 解析
# ============================================================

def parse_pdb(pdb_path: Path) -> dict:
    """解析 PDB 文件"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('complex', str(pdb_path))
    return structure


def get_chain_type(structure) -> dict:
    """识别每个链的类型"""
    chain_info = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            n_residues = len(list(chain.get_residues()))
            
            if n_residues > 150:
                chain_type = 'IL24'
            elif n_residues > 100:
                chain_type = 'VH'
            else:
                chain_type = 'VL'
            
            chain_info[chain_id] = {
                'type': chain_type,
                'n_residues': n_residues
            }
    
    return chain_info


def get_ca_coordinates(structure) -> dict:
    """提取所有 Cα 原子的坐标"""
    coords = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            coords[chain_id] = []
            
            for residue in chain.get_residues():
                if residue.id[0] == ' ':  # 标准氨基酸
                    try:
                        ca = residue['CA']
                        coord = ca.get_coord()
                        coords[chain_id].append({
                            'res_num': residue.id[1],
                            'res_name': PDB.Polypeptide.three_to_one(residue.resname),
                            'coord': coord
                        })
                    except KeyError:
                        continue
    
    return coords


# ============================================================
# 质量评估函数
# ============================================================

def compute_interface_area(structure, chain_info: dict, 
                          distance_threshold: float = 8.0) -> float:
    """
    计算界面面积 (Å²)
    使用简化的 surface accessibility 方法
    """
    # 获取所有原子坐标
    all_atoms = []
    chains = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            chains[chain_id] = []
            
            for residue in chain.get_residues():
                if residue.id[0] == ' ':
                    for atom in residue:
                        all_atoms.append({
                            'chain': chain_id,
                            'res_num': residue.id[1],
                            'coord': atom.get_coord(),
                            'radius': atom.radius if hasattr(atom, 'radius') else 1.4
                        })
                        chains[chain_id].append(all_atoms[-1])
    
    # 计算界面残基
    il24_chains = [c for c, info in chain_info.items() if info['type'] == 'IL24']
    ab_chains = [c for c, info in chain_info.items() if info['type'] in ['VH', 'VL']]
    
    interface_residues = set()
    
    for il24_atom in [a for a in all_atoms if a['chain'] in il24_chains]:
        for ab_atom in [a for a in all_atoms if a['chain'] in ab_chains]:
            dist = np.linalg.norm(il24_atom['coord'] - ab_atom['coord'])
            if dist < distance_threshold:
                interface_residues.add((il24_atom['chain'], il24_atom['res_num']))
                interface_residues.add((ab_atom['chain'], ab_atom['res_num']))
    
    # 简化: 每个残基约 150 Å²
    buried_area = len(interface_residues) * 150
    
    return buried_area, interface_residues


def compute_clashes(structure, chain_info: dict, 
                   clash_threshold: float = 2.0) -> dict:
    """
    检测原子冲突
    
    返回: clash_count, clash_details
    """
    # 获取所有重原子
    heavy_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain.get_residues():
                if residue.id[0] == ' ':
                    for atom in residue:
                        if atom.element != 'H':  # 排除氢原子
                            heavy_atoms.append({
                                'chain': chain.id,
                                'coord': atom.get_coord(),
                                'vdw_radius': get_vdw_radius(atom.element)
                            })
    
    # 检测冲突
    clashes = []
    n_il24 = sum(1 for a in heavy_atoms if chain_info.get(a['chain'], {}).get('type') == 'IL24')
    n_ab = sum(1 for a in heavy_atoms if chain_info.get(a['chain'], {}).get('type') in ['VH', 'VL'])
    
    for i, atom1 in enumerate(heavy_atoms):
        type1 = chain_info.get(atom1['chain'], {}).get('type', '')
        if type1 not in ['IL24', 'VH', 'VL']:
            continue
            
        for j, atom2 in enumerate(heavy_atoms[i+1:], i+1):
            type2 = chain_info.get(atom2['chain'], {}).get('type', '')
            if type2 not in ['IL24', 'VH', 'VL']:
                continue
            
            # 只检查 IL-24 和抗体之间的冲突
            if type1 == type2:
                continue
            
            dist = np.linalg.norm(atom1['coord'] - atom2['coord'])
            min_dist = atom1['vdw_radius'] + atom2['vdw_radius']
            
            if dist < min_dist * 0.8:  # 80% of sum of vdw radii
                clashes.append({
                    'chain1': atom1['chain'],
                    'chain2': atom2['chain'],
                    'distance': dist,
                    'min_expected': min_dist,
                    'overlap': min_dist - dist
                })
    
    return len(clashes), clashes


def get_vdw_radius(element: str) -> float:
    """获取原子的范德华半径 (Å)"""
    vdw_radii = {
        'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
        'H': 1.2, 'P': 1.8, 'FE': 1.35, 'ZN': 1.45,
        'CA': 1.86, 'MG': 1.45, 'MN': 1.45, 'CU': 1.4,
    }
    return vdw_radii.get(element.upper(), 1.5)


def compute_ptm_from_pae(pae_matrix: np.ndarray) -> float:
    """
    从 PAE 矩阵估算 pTM 分数
    这是近似值，实际 pTM 需要从 AF3 获得
    """
    # 简化: 平均 PAE 越低越好
    mean_pae = np.mean(pae_matrix)
    # 假设 PAE 0-30 对应 pTM 1-0
    ptm = max(0, 1 - mean_pae / 30)
    return ptm


def extract_pae_from_pdb(pdb_path: Path) -> np.ndarray:
    """
    从 PDB 文件的 b-factor 列提取 PAE 信息
    AlphaFold 输出中 b-factor 通常是 pLDDT
    """
    pae_data = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                b_factor = float(line[60:66].strip())
                pae_data.append(b_factor)
    
    if not pae_data:
        return np.array([])
    
    return np.array(pae_data)


def compute_interface_contacts(structure, chain_info: dict,
                               distance_threshold: float = 4.5) -> dict:
    """
    计算界面接触数
    """
    il24_chains = [c for c, info in chain_info.items() if info['type'] == 'IL24']
    ab_chains = [c for c, info in chain_info.items() if info['type'] in ['VH', 'VL']]
    
    contacts = []
    contact_types = defaultdict(int)
    
    for model in structure:
        for chain in model:
            if chain.id not in il24_chains:
                continue
            
            for residue1 in chain.get_residues():
                if residue1.id[0] != ' ':
                    continue
                
                for other_chain in model:
                    if other_chain.id not in ab_chains:
                        continue
                    
                    for residue2 in other_chain.get_residues():
                        if residue2.id[0] != ' ':
                            continue
                        
                        # 计算 Cα 距离
                        try:
                            ca1 = residue1['CA'].get_coord()
                            ca2 = residue2['CA'].get_coord()
                            dist = np.linalg.norm(ca1 - ca2)
                            
                            if dist < distance_threshold:
                                contacts.append({
                                    'il24_res': residue1.id[1],
                                    'il24_name': PDB.Polypeptide.three_to_one(residue1.resname),
                                    'ab_chain': other_chain.id,
                                    'ab_res': residue2.id[1],
                                    'ab_name': PDB.Polypeptide.three_to_one(residue2.resname),
                                    'distance': dist
                                })
                                
                                # 统计接触类型
                                contact_types['total'] += 1
                        except KeyError:
                            continue
    
    return contacts, contact_types


def analyze_cdr_involvement(contacts: list, paratope_file: Path = None) -> dict:
    """
    分析哪些 CDR 参与了接触
    """
    # CDR 区域定义 (基于 Kabat)
    cdr_regions = {
        'H': {
            'CDR1': (22, 29),
            'CDR2': (46, 52),
            'CDR3': (95, 105),
        },
        'L': {
            'CDR1': (22, 27),
            'CDR2': (41, 42),
            'CDR3': (85, 90),
        }
    }
    
    involvement = {}
    
    for contact in contacts:
        ab_chain = contact['ab_chain']
        ab_res = contact['ab_res']
        
        chain_type = 'H' if 'H' in ab_chain else 'L'
        
        for cdr_name, (start, end) in cdr_regions.get(chain_type, {}).items():
            if start <= ab_res <= end:
                if cdr_name not in involvement:
                    involvement[cdr_name] = {'count': 0, 'residues': set()}
                involvement[cdr_name]['count'] += 1
                involvement[cdr_name]['residues'].add(ab_res)
    
    return involvement


def compute_epitope_consistency(ensemble_dir: Path) -> dict:
    """
    计算 ensemble 中表位的一致性
    """
    pdb_files = list(ensemble_dir.glob("*.pdb"))
    
    if len(pdb_files) < 2:
        return {'n_models': len(pdb_files), 'consistency': None}
    
    all_epitopes = []
    
    for pdb_file in pdb_files:
        try:
            structure = parse_pdb(pdb_file)
            chain_info = get_chain_type(structure)
            contacts, _ = compute_interface_contacts(structure, chain_info)
            
            epitope_residues = set([c['il24_res'] for c in contacts])
            all_epitopes.append(epitope_residues)
        except Exception as e:
            print(f"⚠️ 处理 {pdb_file.name} 失败: {e}")
    
    if len(all_epitopes) < 2:
        return {'n_models': len(all_epitopes), 'consistency': None}
    
    # 计算每对模型之间的 Jaccard 相似度
    jaccard_scores = []
    for i, ep1 in enumerate(all_epitopes):
        for ep2 in all_epitopes[i+1:]:
            if ep1 or ep2:
                intersection = len(ep1 & ep2)
                union = len(ep1 | ep2)
                jaccard = intersection / union if union > 0 else 0
                jaccard_scores.append(jaccard)
    
    mean_jaccard = np.mean(jaccard_scores) if jaccard_scores else 0
    
    # 计算每个位置出现的频率
    all_residues = set().union(*all_epitopes)
    residue_frequency = {}
    for res in all_residues:
        count = sum(1 for ep in all_epitopes if res in ep)
        residue_frequency[res] = {
            'count': count,
            'frequency': count / len(all_epitopes) * 100
        }
    
    return {
        'n_models': len(all_epitopes),
        'mean_jaccard': mean_jaccard,
        'consistency': 'HIGH' if mean_jaccard > 0.7 else ('MEDIUM' if mean_jaccard > 0.4 else 'LOW'),
        'residue_frequency': residue_frequency
    }


# ============================================================
# 质量报告生成
# ============================================================

def generate_quality_report(metrics: dict, output_dir: Path) -> dict:
    """生成质量报告"""
    
    # 总体评分
    scores = []
    
    # ipTM 评分
    if 'iptm' in metrics:
        iptm = metrics['iptm']
        if iptm >= HIGH_CONFIDENCE_IPTM:
            scores.append(1.0)
        elif iptm >= MEDIUM_CONFIDENCE_IPTM:
            scores.append(0.7)
        else:
            scores.append(0.3)
    
    # 冲突评分
    if 'clash_count' in metrics:
        n_clashes = metrics['clash_count']
        if n_clashes == 0:
            scores.append(1.0)
        elif n_clashes < 10:
            scores.append(0.8)
        elif n_clashes < 50:
            scores.append(0.5)
        else:
            scores.append(0.1)
    
    # 界面面积评分
    if 'interface_area' in metrics:
        area = metrics['interface_area']
        if 1000 <= area <= 3000:
            scores.append(1.0)
        elif 500 <= area < 1000 or 3000 < area <= 5000:
            scores.append(0.7)
        else:
            scores.append(0.4)
    
    # 一致性评分
    if 'consistency' in metrics:
        cons = metrics['consistency']
        if cons == 'HIGH':
            scores.append(1.0)
        elif cons == 'MEDIUM':
            scores.append(0.7)
        else:
            scores.append(0.4)
    
    overall_score = np.mean(scores) if scores else None
    
    # 判定
    if overall_score and overall_score >= 0.8:
        quality = 'EXCELLENT'
    elif overall_score and overall_score >= 0.6:
        quality = 'GOOD'
    elif overall_score and overall_score >= 0.4:
        quality = 'FAIR'
    else:
        quality = 'POOR'
    
    report = {
        'overall_score': overall_score,
        'quality_grade': quality,
        'metrics': metrics,
        'scores': scores,
        'recommendations': []
    }
    
    # 添加建议
    if metrics.get('clash_count', 0) > 10:
        report['recommendations'].append(
            '检测到较多原子冲突，建议进行结构优化'
        )
    
    if metrics.get('interface_area', 0) < 500:
        report['recommendations'].append(
            '界面面积较小，可能表示弱结合'
        )
    
    if metrics.get('consistency') == 'LOW':
        report['recommendations'].append(
            'Ensemble 一致性低，建议增加采样或检查输入'
        )
    
    return report


# ============================================================
# 主函数
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Step 14: Quality Assessment')
    parser.add_argument('--input', '-i', type=str,
                       help='输入目录 (ensemble)',
                       default=None)
    parser.add_argument('--mode', '-m', type=str,
                       choices=['12A_blind', '12B_guided'],
                       default='12A_blind',
                       help='预测模式')
    parser.add_argument('--paratope', '-p', type=str,
                       help='对位文件 (来自 Step 13)',
                       default=None)
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 14: 复合物质量评估")
    print("=" * 60)
    
    # 确定输入路径
    if args.input:
        input_path = Path(args.input)
    else:
        # 尝试自动找到输出目录
        possible_dirs = [
            RESULTS_DIR / "12_complex_prediction" / args.mode,
            RESULTS_DIR / "12_complex_prediction" / args.mode / "af3_output",
        ]
        input_path = None
        for d in possible_dirs:
            if d.exists():
                pdb_files = list(d.glob("*.pdb"))
                if pdb_files:
                    input_path = d
                    break
        
        if input_path is None:
            input_path = RESULTS_DIR / "12_complex_prediction" / args.mode
    
    if not input_path.exists():
        print(f"❌ 错误: 找不到输入 {input_path}")
        sys.exit(1)
    
    print(f"📂 输入: {input_path}")
    
    # 如果是目录，进行 ensemble 分析
    pdb_files = list(input_path.glob("*.pdb"))
    
    if not pdb_files:
        print(f"⚠️ 未找到 PDB 文件")
        sys.exit(1)
    
    print(f"📊 分析 {len(pdb_files)} 个模型...")
    
    all_metrics = []
    
    for pdb_file in pdb_files[:5]:  # 分析前5个模型
        print(f"\n📄 {pdb_file.name}")
        
        try:
            structure = parse_pdb(pdb_file)
            chain_info = get_chain_type(structure)
            
            # 界面面积
            interface_area, interface_res = compute_interface_area(
                structure, chain_info
            )
            print(f"   界面面积: {interface_area:.1f} Å²")
            print(f"   界面残基: {len(interface_res)}")
            
            # 冲突检测
            clash_count, clashes = compute_clashes(structure, chain_info)
            print(f"   原子冲突: {clash_count}")
            
            # 界面接触
            contacts, contact_types = compute_interface_contacts(
                structure, chain_info
            )
            print(f"   界面接触: {len(contacts)}")
            
            # PAE/pLDDT
            pae_data = extract_pae_from_pdb(pdb_file)
            if len(pae_data) > 0:
                mean_plddt = np.mean(pae_data)
                print(f"   平均 pLDDT: {mean_plddt:.1f}")
            
            metrics = {
                'file': pdb_file.name,
                'interface_area': interface_area,
                'n_interface_residues': len(interface_res),
                'clash_count': clash_count,
                'n_contacts': len(contacts),
                'mean_plddt': np.mean(pae_data) if len(pae_data) > 0 else None,
                'chain_info': chain_info
            }
            
            all_metrics.append(metrics)
            
        except Exception as e:
            print(f"   ❌ 分析失败: {e}")
    
    # Ensemble 一致性分析
    print("\n" + "-" * 40)
    print("📊 Ensemble 一致性分析:")
    
    consistency = compute_epitope_consistency(input_path)
    
    if consistency.get('n_models', 0) >= 2:
        print(f"   模型数: {consistency['n_models']}")
        print(f"   平均 Jaccard: {consistency.get('mean_jaccard', 0):.3f}")
        print(f"   一致性: {consistency.get('consistency', 'N/A')}")
    
    # 生成报告
    avg_metrics = {
        'n_models_analyzed': len(all_metrics),
        'mean_interface_area': np.mean([m['interface_area'] for m in all_metrics]),
        'mean_clash_count': np.mean([m['clash_count'] for m in all_metrics]),
        'mean_contacts': np.mean([m['n_contacts'] for m in all_metrics]),
        'ensemble': consistency
    }
    
    report = generate_quality_report(avg_metrics, OUTPUT_DIR)
    
    # 保存报告
    report_file = OUTPUT_DIR / f"quality_report_{args.mode}.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n📋 质量报告: {report_file}")
    
    # 保存指标表
    metrics_df = pd.DataFrame(all_metrics)
    metrics_file = OUTPUT_DIR / f"metrics_{args.mode}.csv"
    metrics_df.to_csv(metrics_file, index=False)
    print(f"📊 指标表: {metrics_file}")
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("📊 质量评估摘要")
    print("=" * 60)
    print(f"   总体评分: {report['overall_score']:.2f}" if report['overall_score'] else "   N/A")
    print(f"   质量等级: {report['quality_grade']}")
    print(f"   平均界面面积: {avg_metrics['mean_interface_area']:.1f} Å²")
    print(f"   平均冲突数: {avg_metrics['mean_clash_count']:.1f}")
    print(f"   平均接触数: {avg_metrics['mean_contacts']:.1f}")
    
    if report['recommendations']:
        print("\n📋 建议:")
        for rec in report['recommendations']:
            print(f"   - {rec}")
    
    print("\n" + "=" * 60)
    print("Step 14 完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
