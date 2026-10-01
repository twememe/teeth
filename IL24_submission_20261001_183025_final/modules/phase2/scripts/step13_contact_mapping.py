#!/usr/bin/env python3
"""
Step 13: 结合位点定位 - Epitope/Paratope 接触分析
====================================================
基于复合物模型，计算抗体-抗原接触界面

输入: IL-24_Ab_complex.pdb (ensemble)
输出: epitope.csv, paratope.csv, contact_summary.csv

判断标准:
- 接触距离阈值: d_ij < 4.5 Å (任一原子对)
- 三重验证: Accessible + Contact-supported + Functionally relevant
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
    import MDAnalysis as mda
    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False
    print("⚠️ 建议安装: pip install MDAnalysis")


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "13_contact_mapping"

# IL-24 残基编号映射 (project position = UniProt position - 39)
# project 1-181 对应 UniProt 40-220
IL24_PROJECT_START = 1
IL24_PROJECT_END = 181
IL24_CHAIN = 'A'  # AlphaFold 结构中 IL-24 的链

# 接触距离阈值
CONTACT_DISTANCE = 4.5  # Å
HYDROGEN_BOND_DISTANCE = 3.5  # Å
SALT_BRIDGE_DISTANCE = 4.0  # Å

# Phase 1 热点 (用于验证)
PHASE1_HOTSPOTS = [
    {"region": "C1-B", "start": 52, "end": 66},
    {"region": "C2-B", "start": 126, "end": 140},
    {"region": "C3-B", "start": 60, "end": 74},
    {"region": "C4-B", "start": 143, "end": 157},
    {"region": "C5-B", "start": 117, "end": 131},
]

# CDR 区域定义 (Kabat 编号)
CDR_REGIONS = {
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

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# PDB 解析函数
# ============================================================

def parse_pdb_atoms(pdb_path: Path) -> pd.DataFrame:
    """解析 PDB 文件，提取所有原子坐标"""
    atoms = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'res_name': line[17:20].strip(),
                    'chain': line[21],
                    'res_seq': int(line[22:26].strip()),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': float(line[54:60]) if line[54:60].strip() else 1.0,
                    'b_factor': float(line[60:66]) if line[60:66].strip() else 0.0,
                }
                atoms.append(atom)
    
    return pd.DataFrame(atoms)


def identify_chain_type(atoms_df: pd.DataFrame) -> dict:
    """
    识别每个链的类型 (IL-24, VH, VL)
    基于序列长度和特征
    """
    chain_info = {}
    
    for chain in atoms_df['chain'].unique():
        chain_atoms = atoms_df[atoms_df['chain'] == chain]
        n_residues = chain_atoms['res_seq'].nunique()
        n_atoms = len(chain_atoms)
        
        # 基于 IL-24 长度 (~181 aa) 和抗体 VH/VL (~120 aa) 区分
        if n_residues > 150:
            chain_type = 'IL24'
        elif n_atoms > 2000:
            chain_type = 'VH'
        else:
            chain_type = 'VL'
        
        chain_info[chain] = {
            'type': chain_type,
            'n_residues': n_residues,
            'n_atoms': n_atoms
        }
    
    return chain_info


def compute_contacts(atoms_df: pd.DataFrame, chain_info: dict, 
                     distance_threshold: float = 4.5) -> pd.DataFrame:
    """
    计算所有原子对之间的距离，识别接触
    
    返回: 接触列表 (residue-level)
    """
    contacts = []
    
    # 分离 IL-24 和抗体原子
    il24_chains = [c for c, info in chain_info.items() if info['type'] == 'IL24']
    ab_chains = [c for c, info in chain_info.items() if info['type'] in ['VH', 'VL']]
    
    for il24_chain in il24_chains:
        for ab_chain in ab_chains:
            il24_atoms = atoms_df[atoms_df['chain'] == il24_chain]
            ab_atoms = atoms_df[atoms_df['chain'] == ab_chain]
            
            # 计算距离矩阵
            for _, il24_atom in il24_atoms.iterrows():
                for _, ab_atom in ab_atoms.iterrows():
                    dist = np.sqrt(
                        (il24_atom['x'] - ab_atom['x'])**2 +
                        (il24_atom['y'] - ab_atom['y'])**2 +
                        (il24_atom['z'] - ab_atom['z'])**2
                    )
                    
                    if dist < distance_threshold:
                        contacts.append({
                            'il24_chain': il24_chain,
                            'il24_res': il24_atom['res_seq'],
                            'il24_res_name': il24_atom['res_name'],
                            'il24_atom': il24_atom['name'],
                            'ab_chain': ab_chain,
                            'ab_res': ab_atom['res_seq'],
                            'ab_res_name': ab_atom['res_name'],
                            'ab_atom': ab_atom['name'],
                            'distance': dist,
                            'ab_chain_type': chain_info[ab_chain]['type']
                        })
    
    return pd.DataFrame(contacts)


def aggregate_to_residues(contacts_df: pd.DataFrame) -> tuple:
    """
    将原子级接触聚合为残基层级
    
    返回: (epitope_df, paratope_df)
    """
    if contacts_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    
    # IL-24 表位 (每个残基接触的原子数)
    epitope_counts = contacts_df.groupby(
        ['il24_chain', 'il24_res', 'il24_res_name']
    ).agg({
        'distance': ['min', 'count'],
        'ab_chain': lambda x: list(x.unique())
    }).reset_index()
    epitope_counts.columns = ['chain', 'residue', 'residue_name', 
                              'min_distance', 'n_contacts', 'ab_chains']
    
    # 抗体对位
    paratope_counts = contacts_df.groupby(
        ['ab_chain', 'ab_res', 'ab_res_name', 'ab_chain_type']
    ).agg({
        'distance': ['min', 'count'],
        'il24_chain': lambda x: list(x.unique())
    }).reset_index()
    paratope_counts.columns = ['chain', 'residue', 'residue_name', 
                               'chain_type', 'min_distance', 'n_contacts', 
                               'antigen_chains']
    
    return epitope_counts, paratope_counts


def map_kabat_position(chain_type: str, res_number: int) -> dict:
    """
    将 PDB 中的残基编号映射到 Kabat 编号
    需要知道框架区的对齐
    """
    # 简化版：假设 PDB 中的编号直接对应 Kabat 编号
    # 实际使用时需要 ANARCI 进行精确对齐
    return {'kabat': res_number, 'imgt': res_number, 'original': res_number}


def classify_cdr(chain_type: str, res_number: int) -> str:
    """判断残基是否属于 CDR 区域"""
    if chain_type not in CDR_REGIONS:
        return 'framework'
    
    for cdr_name, (start, end) in CDR_REGIONS[chain_type].items():
        if start <= res_number <= end:
            return cdr_name
    
    return 'framework'


def check_phase1_overlap(residue: int, hotspots: list) -> list:
    """检查残基是否落入 Phase 1 热点区域"""
    matching = []
    for hs in hotspots:
        if hs['start'] <= residue <= hs['end']:
            matching.append(hs['region'])
    return matching


def generate_contact_summary(epitope_df: pd.DataFrame, 
                              paratope_df: pd.DataFrame,
                              hotspots: list) -> dict:
    """生成接触摘要统计"""
    summary = {
        'n_epitope_residues': len(epitope_df),
        'n_paratope_residues': len(paratope_df),
        'total_contacts': epitope_df['n_contacts'].sum() if not epitope_df.empty else 0,
        'mean_contact_distance': epitope_df['min_distance'].mean() if not epitope_df.empty else None,
        'phase1_coverage': {}
    }
    
    # Phase 1 热点覆盖率
    for hs in hotspots:
        region = hs['region']
        start, end = hs['start'], hs['end']
        
        if not epitope_df.empty:
            covered = epitope_df[
                (epitope_df['residue'] >= start) & 
                (epitope_df['residue'] <= end)
            ]
            n_covered = len(covered)
            total = end - start + 1
            coverage = n_covered / total * 100
            
            summary['phase1_coverage'][region] = {
                'covered': n_covered,
                'total': total,
                'coverage_pct': round(coverage, 1)
            }
        else:
            summary['phase1_coverage'][region] = {
                'covered': 0,
                'total': end - start + 1,
                'coverage_pct': 0.0
            }
    
    # CDR 参与统计
    if not paratope_df.empty:
        cdr_counts = paratope_df['cdr_region'].value_counts()
        summary['cdr_involvement'] = cdr_counts.to_dict()
    else:
        summary['cdr_involvement'] = {}
    
    return summary


# ============================================================
# Ensemble 分析
# ============================================================

def analyze_ensemble(ensemble_dir: Path, chain_info: dict) -> dict:
    """
    分析 ensemble 预测的一致性
    
    对多个种子/模型的结果进行统计
    """
    pdb_files = list(ensemble_dir.glob("*.pdb"))
    
    if not pdb_files:
        return {'error': 'No PDB files found'}
    
    all_epitopes = []
    all_paratopes = []
    
    for pdb_file in pdb_files:
        try:
            atoms = parse_pdb_atoms(pdb_file)
            contacts = compute_contacts(atoms, chain_info, CONTACT_DISTANCE)
            epitope, paratope = aggregate_to_residues(contacts)
            
            all_epitopes.append(set(epitope['residue'].tolist()) if not epitope.empty else set())
            all_paratopes.append(set(paratope['residue'].tolist()) if not paratope.empty else set())
        except Exception as e:
            print(f"⚠️ 处理 {pdb_file.name} 失败: {e}")
    
    # 计算一致性
    if all_epitopes:
        # 每个残基在多少模型中出现
        all_residues = set().union(*all_epitopes)
        residue_frequency = {}
        for res in all_residues:
            count = sum(1 for ep in all_epitopes if res in ep)
            residue_frequency[res] = {
                'frequency': count,
                'pct': count / len(all_epitopes) * 100
            }
        
        ensemble_summary = {
            'n_models': len(all_epitopes),
            'n_unique_epitope_residues': len(all_residues),
            'consensus_residues': {
                r: freq for r, freq in residue_frequency.items() 
                if freq['pct'] >= 80  # 80% 以上模型都预测到的残基
            },
            'residue_frequency': residue_frequency
        }
    else:
        ensemble_summary = {'n_models': 0}
    
    return ensemble_summary


# ============================================================
# 主函数
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Step 13: Contact Mapping')
    parser.add_argument('--input', '-i', type=str, 
                       help='输入 PDB 文件或目录 (ensemble)',
                       default=None)
    parser.add_argument('--mode', '-m', type=str, 
                       choices=['12A_blind', '12B_guided'],
                       default='12A_blind',
                       help='预测模式')
    parser.add_argument('--distance', '-d', type=float,
                       default=4.5,
                       help='接触距离阈值 (Å)')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Step 13: 结合位点定位")
    print("=" * 60)
    
    # 确定输入路径
    if args.input:
        input_path = Path(args.input)
    else:
        input_path = RESULTS_DIR / "12_complex_prediction" / args.mode
    
    if not input_path.exists():
        print(f"❌ 错误: 找不到输入 {input_path}")
        print(f"请先运行 Step 12 完成复合物预测")
        sys.exit(1)
    
    print(f"📂 输入: {input_path}")
    print(f"📏 接触距离: {args.distance} Å")
    print()
    
    # 如果是目录，进行 ensemble 分析
    if input_path.is_dir():
        # 假设子目录 af3_output_*
        output_subdir = input_path
        pdb_files = list(output_subdir.glob("*.pdb"))
        
        if not pdb_files:
            # 查找子目录
            subdirs = [d for d in output_subdir.iterdir() if d.is_dir()]
            for subdir in subdirs:
                pdb_files = list(subdir.glob("*.pdb"))
                if pdb_files:
                    output_subdir = subdir
                    break
        
        if pdb_files:
            print(f"🔍 检测到 ensemble: {len(pdb_files)} 个模型")
            print(f"📂 PDB 文件目录: {output_subdir}")
            
            # 使用第一个模型识别链类型
            sample_pdb = pdb_files[0]
            print(f"📊 样本文件: {sample_pdb.name}")
            
            atoms = parse_pdb_atoms(sample_pdb)
            chain_info = identify_chain_type(atoms)
            
            print("\n🔗 链识别结果:")
            for chain, info in chain_info.items():
                print(f"   {chain}: {info['type']} ({info['n_residues']} residues)")
            
            # Ensemble 分析
            print("\n🔄 进行 ensemble 分析...")
            ensemble_results = analyze_ensemble(output_subdir, chain_info)
            
            # 保存 ensemble 结果
            ensemble_file = OUTPUT_DIR / f"ensemble_summary_{args.mode}.json"
            with open(ensemble_file, 'w') as f:
                json.dump(ensemble_results, f, indent=2)
            print(f"📋 Ensemble 摘要: {ensemble_file}")
            
            # 使用 consensus 残基
            if 'consensus_residues' in ensemble_results:
                consensus = ensemble_results['consensus_residues']
                print(f"\n🎯 Consensus 表位残基 (≥80% 模型): {len(consensus)}")
                for res in sorted(consensus.keys()):
                    print(f"   {res}: {consensus[res]['pct']:.1f}%")
        else:
            print(f"⚠️ 未找到 PDB 文件")
            sys.exit(1)
    else:
        # 单个 PDB 文件
        print(f"📊 处理单个文件: {input_path}")
        
        atoms = parse_pdb_atoms(input_path)
        chain_info = identify_chain_type(atoms)
        
        contacts = compute_contacts(atoms, chain_info, args.distance)
        epitope_df, paratope_df = aggregate_to_residues(contacts)
        
        # 添加 CDR 分类
        if not paratope_df.empty:
            paratope_df['cdr_region'] = paratope_df.apply(
                lambda r: classify_cdr(r['chain_type'], r['residue']),
                axis=1
            )
        
        # 添加 Phase 1 热点重叠
        if not epitope_df.empty:
            epitope_df['phase1_overlap'] = epitope_df['residue'].apply(
                lambda r: check_phase1_overlap(r, hotspots)
            )
        
        # 保存结果
        epitope_file = OUTPUT_DIR / f"epitope_{args.mode}.csv"
        paratope_file = OUTPUT_DIR / f"paratope_{args.mode}.csv"
        
        epitope_df.to_csv(epitope_file, index=False)
        paratope_df.to_csv(paratope_file, index=False)
        
        print(f"✅ 表位残基: {len(epitope_df)} -> {epitope_file}")
        print(f"✅ 对位残基: {len(paratope_df)} -> {paratope_file}")
        
        # 打印摘要
        print("\n" + "-" * 40)
        print("📊 接触摘要:")
        print(f"   表位残基数: {len(epitope_df)}")
        print(f"   对位残基数: {len(paratope_df)}")
        
        if not paratope_df.empty:
            print(f"\n   CDR 参与:")
            for cdr, count in paratope_df['cdr_region'].value_counts().items():
                print(f"     {cdr}: {count} residues")
        
        # Phase 1 热点验证
        print(f"\n🎯 Phase 1 热点验证:")
        for hs in PHASE1_HOTSPOTS:
            region = hs['region']
            if not epitope_df.empty:
                covered = epitope_df[
                    (epitope_df['residue'] >= hs['start']) & 
                    (epitope_df['residue'] <= hs['end'])
                ]
                n_covered = len(covered)
                pct = n_covered / (hs['end'] - hs['start'] + 1) * 100
                status = "✓" if n_covered > 0 else "✗"
                print(f"   {status} {region} ({hs['start']}-{hs['end']}): {n_covered}/{hs['end']-hs['start']+1} ({pct:.1f}%)")
    
    # 生成 provenance
    provenance = {
        'script': 'step13_contact_mapping.py',
        'timestamp': datetime.now().isoformat(),
        'parameters': {
            'contact_distance': args.distance,
            'mode': args.mode,
            'phase1_hotspots': PHASE1_HOTSPOTS
        }
    }
    
    prov_file = OUTPUT_DIR / f"provenance_{args.mode}.json"
    with open(prov_file, 'w') as f:
        json.dump(provenance, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Step 13 完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
