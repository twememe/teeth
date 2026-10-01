#!/usr/bin/env python3
"""
Step 16: 竞争性结合分析
=========================
判断抗体中和机制：直接阻断 vs 变构干扰

分析内容:
- E_Ab vs E_Receptor 重叠度
- Jaccard 指数
- Receptor-contact coverage
- Steric clash 检测

输入: IL-24_Ab_complex.pdb, receptor_binding_sites.json
输出: competition_report.json
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
OUTPUT_DIR = RESULTS_DIR / "16_competition_analysis"

# IL-24 受体结合位点 (基于文献)
# IL-24 主要与 IL-20R1/IL-20R2 和 IL-22R1/IL-20R2 受体复合物结合
RECEPTOR_BINDING_SITES = {
    # IL-20R1 结合区域
    'IL20R1_region1': list(range(55, 75)),
    'IL20R1_region2': list(range(120, 140)),
    # IL-22R1 结合区域
    'IL22R1_region1': list(range(65, 85)),
    'IL22R1_region2': list(range(130, 150)),
    # 通用受体结合面
    'core_interface': list(range(60, 80)) + list(range(125, 145)),
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 重叠度计算
# ============================================================

def compute_overlap(set_a: set, set_b: set) -> dict:
    """
    计算两个集合的重叠
    
    返回:
    - intersection: 交集
    - union: 并集
    - overlap_count: 交集大小
    - jaccard: Jaccard 指数
    - overlap_coefficient: 重叠系数
    """
    intersection = set_a & set_b
    union = set_a | set_b
    
    overlap_count = len(intersection)
    jaccard = len(intersection) / len(union) if len(union) > 0 else 0
    overlap_coefficient = (
        len(intersection) / min(len(set_a), len(set_b)) 
        if min(len(set_a), len(set_b)) > 0 else 0
    )
    
    return {
        'intersection': sorted(list(intersection)),
        'union': sorted(list(union)),
        'overlap_count': overlap_count,
        'set_a_size': len(set_a),
        'set_b_size': len(set_b),
        'jaccard': round(jaccard, 4),
        'overlap_coefficient': round(overlap_coefficient, 4)
    }


def compute_receptor_contact_coverage(antibody_epitope: set,
                                     receptor_sites: dict) -> dict:
    """
    计算抗体对受体接触面的覆盖程度
    
    返回每个受体区域的覆盖情况
    """
    coverage = {}
    
    for region_name, residues in receptor_sites.items():
        receptor_set = set(residues)
        overlap = antibody_epitope & receptor_set
        coverage_rate = len(overlap) / len(receptor_set) if len(receptor_set) > 0 else 0
        
        coverage[region_name] = {
            'receptor_residues': len(receptor_set),
            'covered_by_antibody': len(overlap),
            'coverage_rate': round(coverage_rate, 4),
            'covered_positions': sorted(list(overlap))
        }
    
    # 计算总体覆盖
    all_receptor = set()
    for residues in receptor_sites.values():
        all_receptor.update(residues)
    
    total_overlap = antibody_epitope & all_receptor
    total_coverage = len(total_overlap) / len(all_receptor) if len(all_receptor) > 0 else 0
    
    coverage['overall'] = {
        'total_receptor_residues': len(all_receptor),
        'covered_by_antibody': len(total_overlap),
        'coverage_rate': round(total_coverage, 4)
    }
    
    return coverage


def detect_steric_clash(antibody_epitope: set, 
                       receptor_sites: dict,
                       interface_area: float) -> dict:
    """
    检测空间位阻
    
    如果抗体完全覆盖受体结合面，则存在空间位阻风险
    """
    # 合并所有受体位点
    all_receptor = set()
    for residues in receptor_sites.values():
        all_receptor.update(residues)
    
    # 计算重叠
    overlap = antibody_epitope & all_receptor
    overlap_ratio = len(overlap) / len(all_receptor) if len(all_receptor) > 0 else 0
    
    # 判断
    # 如果界面面积大（>2500 Å²）且高度重叠，提示可能存在位阻
    has_clash_risk = (
        overlap_ratio > 0.5 and 
        interface_area > 2500
    )
    
    return {
        'overlap_with_receptor': len(overlap),
        'overlap_ratio': round(overlap_ratio, 4),
        'interface_area': interface_area,
        'steric_clash_risk': has_clash_risk,
        'interpretation': (
            'HIGH_RISK' if overlap_ratio > 0.7 else
            'MODERATE_RISK' if overlap_ratio > 0.5 else
            'LOW_RISK'
        )
    }


def determine_mechanism(overlap_data: dict, 
                       coverage_data: dict,
                       clash_data: dict) -> dict:
    """
    综合判断中和机制
    
    机制分类:
    1. Competitive blocking: E_Ab 与 E_Receptor 高度重叠
    2. Steric hindrance: E_Ab 部分重叠，但覆盖受体路径
    3. Allosteric interference: E_Ab 与 E_Receptor 无重叠
    """
    jaccard = overlap_data['jaccard']
    coverage_rate = coverage_data['overall']['coverage_rate']
    clash_risk = clash_data['steric_clash_risk']
    
    # 判定规则（连续谱，非硬阈值）
    mechanism_scores = {
        'competitive_blocking': 0,
        'steric_hindrance': 0,
        'allosteric_interference': 0
    }
    
    # Jaccard 评分
    if jaccard > 0.5:
        mechanism_scores['competitive_blocking'] += jaccard * 2
        mechanism_scores['steric_hindrance'] += jaccard
    elif jaccard > 0.2:
        mechanism_scores['competitive_blocking'] += jaccard
        mechanism_scores['steric_hindrance'] += jaccard * 2
    else:
        mechanism_scores['allosteric_interference'] += (1 - jaccard) * 2
    
    # 覆盖率评分
    if coverage_rate > 0.5:
        mechanism_scores['competitive_blocking'] += coverage_rate
        mechanism_scores['steric_hindrance'] += coverage_rate * 0.5
    elif coverage_rate > 0.2:
        mechanism_scores['steric_hindrance'] += coverage_rate * 2
    else:
        mechanism_scores['allosteric_interference'] += (1 - coverage_rate)
    
    # 位阻风险评分
    if clash_risk:
        mechanism_scores['steric_hindrance'] += 0.5
    
    # 归一化
    total = sum(mechanism_scores.values())
    if total > 0:
        normalized_scores = {k: v/total for k, v in mechanism_scores.items()}
    else:
        normalized_scores = mechanism_scores
    
    # 判定主要机制
    primary_mechanism = max(normalized_scores, key=normalized_scores.get)
    
    return {
        'mechanism_scores': {k: round(v, 4) for k, v in normalized_scores.items()},
        'primary_mechanism': primary_mechanism,
        'confidence': 'HIGH' if max(normalized_scores.values()) > 0.6 else (
            'MEDIUM' if max(normalized_scores.values()) > 0.4 else 'LOW'
        ),
        'details': {
            'jaccard': jaccard,
            'coverage_rate': coverage_rate,
            'clash_risk': clash_risk
        }
    }


def read_epitope_from_csv(epitope_file: Path) -> set:
    """从 Step 13 的输出读取表位"""
    if not epitope_file.exists():
        return set()
    
    df = pd.read_csv(epitope_file)
    return set(df['residue'].tolist())


def generate_epitope_residue_positions(pdb_path: Path) -> set:
    """从 PDB 文件提取 IL-24 表位残基位置"""
    residues = set()
    
    with open(pdb_path, 'r') as f:
        current_chain = None
        current_res = None
        
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                res = int(line[22:26].strip())
                
                # 假设 IL-24 在链 A
                if chain == 'A':
                    if res != current_res:
                        residues.add(res)
                        current_res = res
    
    return residues


# ============================================================
# 主函数
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Step 16: Competition Analysis')
    parser.add_argument('--epitope', '-e', type=str,
                       help='表位 CSV 文件 (来自 Step 13)',
                       default=None)
    parser.add_argument('--pdb', '-p', type=str,
                       help='复合物 PDB 文件',
                       default=None)
    parser.add_argument('--mode', '-m', type=str,
                       choices=['12A_blind', '12B_guided'],
                       default='12A_blind')
    parser.add_argument('--interface-area', '-a', type=float,
                       help='界面面积 (Å²)',
                       default=2000)
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 16: 竞争性结合分析")
    print("=" * 60)
    
    # 获取抗体表位
    if args.epitope and Path(args.epitope).exists():
        antibody_epitope = read_epitope_from_csv(Path(args.epitope))
        print(f"📂 从 CSV 读取表位: {len(antibody_epitope)} 残基")
    elif args.pdb and Path(args.pdb).exists():
        antibody_epitope = generate_epitope_residue_positions(Path(args.pdb))
        print(f"📂 从 PDB 提取表位: {len(antibody_epitope)} 残基")
    else:
        # 使用 Phase 1 热点作为默认值
        antibody_epitope = set()
        for hs_start, hs_end in [(52, 66), (60, 74), (126, 140), (143, 157), (117, 131)]:
            antibody_epitope.update(range(hs_start, hs_end + 1))
        print(f"📂 使用 Phase 1 热点: {len(antibody_epitope)} 残基")
    
    print(f"\n🎯 抗体表位: {sorted(antibody_epitope)}")
    
    # 获取受体结合位点
    print("\n🧬 受体结合位点:")
    for region, residues in RECEPTOR_BINDING_SITES.items():
        print(f"   {region}: {residues}")
    
    # 合并受体位点
    receptor_sites = RECEPTOR_BINDING_SITES
    
    # Step 1: 计算重叠
    print("\n" + "-" * 40)
    print("📊 重叠度分析")
    
    all_receptor = set()
    for residues in receptor_sites.values():
        all_receptor.update(residues)
    
    overlap_data = compute_overlap(antibody_epitope, all_receptor)
    
    print(f"   抗体表位大小: {overlap_data['set_a_size']}")
    print(f"   受体位点大小: {overlap_data['set_b_size']}")
    print(f"   重叠残基: {overlap_data['overlap_count']}")
    print(f"   Jaccard 指数: {overlap_data['jaccard']}")
    print(f"   重叠系数: {overlap_data['overlap_coefficient']}")
    print(f"   重叠位置: {overlap_data['intersection']}")
    
    # Step 2: 计算受体接触覆盖率
    print("\n" + "-" * 40)
    print("📊 受体接触覆盖率")
    
    coverage_data = compute_receptor_contact_coverage(antibody_epitope, receptor_sites)
    
    for region, data in coverage_data.items():
        if region != 'overall':
            print(f"   {region}: {data['covered_by_antibody']}/{data['receptor_residues']} "
                  f"({data['coverage_rate']*100:.1f}%)")
    
    print(f"   ─────────────────")
    print(f"   总体覆盖: {coverage_data['overall']['covered_by_antibody']}/"
          f"{coverage_data['overall']['total_receptor_residues']} "
          f"({coverage_data['overall']['coverage_rate']*100:.1f}%)")
    
    # Step 3: 位阻检测
    print("\n" + "-" * 40)
    print("📊 空间位阻分析")
    
    clash_data = detect_steric_clash(antibody_epitope, receptor_sites, args.interface_area)
    
    print(f"   与受体重叠: {clash_data['overlap_with_receptor']}")
    print(f"   重叠比例: {clash_data['overlap_ratio']*100:.1f}%")
    print(f"   界面面积: {clash_data['interface_area']:.1f} Å²")
    print(f"   位阻风险: {clash_data['interpretation']}")
    
    # Step 4: 机制判定
    print("\n" + "-" * 40)
    print("🎯 中和机制判定")
    
    mechanism_data = determine_mechanism(overlap_data, coverage_data, clash_data)
    
    print(f"   主要机制: {mechanism_data['primary_mechanism'].upper()}")
    print(f"   置信度: {mechanism_data['confidence']}")
    print(f"\n   机制评分:")
    for mech, score in mechanism_data['mechanism_scores'].items():
        bar = "█" * int(score * 20)
        print(f"      {mech}: {score:.3f} {bar}")
    
    # 生成详细报告
    report = {
        'script': 'step16_competition_analysis.py',
        'timestamp': datetime.now().isoformat(),
        'mode': args.mode,
        'antibody_epitope': sorted(list(antibody_epitope)),
        'receptor_binding_sites': RECEPTOR_BINDING_SITES,
        'overlap_analysis': overlap_data,
        'receptor_coverage': coverage_data,
        'steric_clash': clash_data,
        'mechanism': mechanism_data,
        'biological_interpretation': {
            'competitive_blocking': (
                '抗体直接阻断 IL-24 与受体结合，'
                '表位与受体结合面高度重叠'
            ),
            'steric_hindrance': (
                '抗体覆盖受体接近 IL-24 的空间路径，'
                '部分重叠，可能通过位阻效应干扰结合'
            ),
            'allosteric_interference': (
                '抗体结合位置远离受体结合面，'
                '可能通过构象变化间接影响受体结合'
            )
        }
    }
    
    # 保存报告
    report_file = OUTPUT_DIR / f"competition_report_{args.mode}.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n📋 报告已保存: {report_file}")
    
    # 打印最终结论
    print("\n" + "=" * 60)
    print("📊 最终结论")
    print("=" * 60)
    
    mech = mechanism_data['primary_mechanism']
    conf = mechanism_data['confidence']
    
    if mech == 'competitive_blocking':
        conclusion = (
            "预测机制: 直接竞争性阻断\n"
            f"置信度: {conf}\n"
            "解读: 抗体表位与 IL-24 受体结合面高度重叠，"
            "通过物理阻断实现中和"
        )
    elif mech == 'steric_hindrance':
        conclusion = (
            "预测机制: 空间位阻\n"
            f"置信度: {conf}\n"
            "解读: 抗体覆盖受体结合面的部分区域，"
            "可能通过位阻效应干扰受体接近"
        )
    else:
        conclusion = (
            "预测机制: 变构干扰\n"
            f"置信度: {conf}\n"
            "解读: 抗体结合位置与受体结合面分离，"
            "可能通过远程构象变化影响功能"
        )
    
    print(conclusion)
    
    print("\n⚠️ 注意: 以上均为计算预测，最终机制需湿实验验证")
    
    print("\n" + "=" * 60)
    print("Step 16 完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
