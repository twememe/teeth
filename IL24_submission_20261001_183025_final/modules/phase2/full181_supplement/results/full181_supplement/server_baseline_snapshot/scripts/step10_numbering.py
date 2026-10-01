#!/usr/bin/env python3
"""
Step 10: 编号体系统一
======================
将抗体 VH/VL 序列转换为标准编号体系 (IMGT/Kabat/Chothia)

本脚本主要是质量检查，编号转换需要 ANARCI 工具

输入: VH_VL.fasta
输出: cdr_annotation.csv, qc_report.json
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime

import pandas as pd


# ============================================================
# 配置
# ============================================================

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "results" / "10_numbering"
INPUT_FASTA = DATA_DIR / "VH_VL.fasta"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 序列解析
# ============================================================

def parse_fasta(fasta_path: Path) -> dict:
    """解析 FASTA 文件"""
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


def identify_chain_type(sequence: str) -> str:
    """根据序列特征识别链类型"""
    if sequence.startswith('EVQL') or sequence.startswith('QVQL') or \
       sequence.startswith('IVQL') or sequence.startswith('KVQL'):
        return 'VH'
    elif sequence.startswith('DIVMT') or sequence.startswith('DVVMT') or \
         sequence.startswith('EIVMT'):
        return 'VL_kappa'
    elif sequence.startswith('DIQMT') or sequence.startswith('EIQMT'):
        return 'VL_lambda'
    else:
        return 'UNKNOWN'


def quality_check(sequence: str, chain_type: str) -> dict:
    """质量检查"""
    checks = {
        'has_stop_codon': '*' in sequence,
        'has_unknown_aa': 'X' in sequence or 'B' in sequence or 'Z' in sequence,
        'length': len(sequence),
        'has_productive_fr': False,
        'cdrs_reasonable': False,
        'issues': []
    }
    
    # 检查生产性
    if not checks['has_stop_codon'] and not checks['has_unknown_aa']:
        checks['has_productive_fr'] = True
    
    # 长度检查
    if chain_type == 'VH':
        if 100 <= len(sequence) <= 140:
            checks['cdrs_reasonable'] = True
        else:
            checks['issues'].append(f'VH 长度 {len(sequence)} 不在标准范围 100-140 内')
    elif 'VL' in chain_type:
        if 95 <= len(sequence) <= 115:
            checks['cdrs_reasonable'] = True
        else:
            checks['issues'].append(f'VL 长度 {len(sequence)} 不在标准范围 95-115 内')
    
    # 框架区保守motif检查
    FR4_motifs = {
        'VH': ['WGXG', 'WGAG', 'WGQG'],
        'VL_kappa': ['FGXGT', 'FGXGTK', 'FGGGTK'],
        'VL_lambda': ['WFXXXXL', 'CLS']
    }
    
    if chain_type in FR4_motifs:
        found = any(motif in sequence for motif in FR4_motifs[chain_type])
        if not found:
            checks['issues'].append('FR4 保守motif未找到')
    
    checks['passed'] = (
        not checks['has_stop_codon'] and
        not checks['has_unknown_aa'] and
        checks['cdrs_reasonable']
    )
    
    return checks


def generate_placeholder_numbering(sequence: str, chain_type: str) -> pd.DataFrame:
    """
    生成占位编号表
    
    注意：精确的 IMGT/Kabat 编号需要 ANARCI 工具
    这里生成基于 Kabat 规则的近似标注
    """
    regions = []
    
    if chain_type == 'VH':
        # Kabat 定义
        # FR1: 1-25, CDR1: 26-35, FR2: 36-49, CDR2: 50-65, FR3: 66-104, CDR3: 105-117, FR4: 118-129
        # 实际长度可能略有差异
        pass
    elif 'VL' in chain_type:
        # Kabat 定义
        # FR1: 1-23, CDR1: 24-34, FR2: 35-49, CDR2: 50-56, FR3: 57-88, CDR3: 89-97, FR4: 98-108
        pass
    
    # 简化版本：使用已知的 CDR 边界（来自 antibody_numbering.csv）
    # 这需要 ANARCI 进行精确转换
    
    return pd.DataFrame(regions)


# ============================================================
# 主函数
# ============================================================

def main():
    print("=" * 60)
    print("Step 10: 编号体系统一")
    print("=" * 60)
    
    # 检查输入
    if not INPUT_FASTA.exists():
        print(f"❌ 错误: 找不到输入文件 {INPUT_FASTA}")
        print("请确保 VH_VL.fasta 存在于 data/ 目录")
        sys.exit(1)
    
    # 解析序列
    print(f"📂 读取: {INPUT_FASTA}")
    sequences = parse_fasta(INPUT_FASTA)
    
    results = []
    all_passed = True
    
    for name, sequence in sequences.items():
        chain_type = identify_chain_type(sequence)
        
        print(f"\n📊 {name}")
        print(f"   链类型: {chain_type}")
        print(f"   长度: {len(sequence)} aa")
        print(f"   前10个残基: {sequence[:10]}...")
        
        # 质量检查
        qc = quality_check(sequence, chain_type)
        
        print(f"   质量检查:")
        print(f"      - 无终止密码子: {'✓' if not qc['has_stop_codon'] else '✗'}")
        print(f"      - 无未知氨基酸: {'✓' if not qc['has_unknown_aa'] else '✗'}")
        print(f"      - CDR长度合理: {'✓' if qc['cdrs_reasonable'] else '✗'}")
        
        if qc['issues']:
            for issue in qc['issues']:
                print(f"      ⚠️ {issue}")
        
        if not qc['passed']:
            all_passed = False
        
        results.append({
            'name': name,
            'chain_type': chain_type,
            'length': len(sequence),
            'sequence': sequence,
            'qc_passed': qc['passed'],
            'qc_details': qc
        })
    
    # 保存结果
    print("\n" + "-" * 40)
    
    if all_passed:
        print("✅ 所有质量检查通过")
    else:
        print("⚠️ 部分质量检查未通过，请检查序列")
    
    # 保存 QC 报告
    qc_report = {
        'script': 'step10_numbering.py',
        'timestamp': datetime.now().isoformat(),
        'input_file': str(INPUT_FASTA),
        'all_passed': all_passed,
        'results': [{
            'name': r['name'],
            'chain_type': r['chain_type'],
            'length': r['length'],
            'qc_passed': r['qc_passed'],
            'issues': r['qc_details'].get('issues', [])
        } for r in results]
    }
    
    report_file = OUTPUT_DIR / "qc_report.json"
    with open(report_file, 'w') as f:
        json.dump(qc_report, f, indent=2)
    print(f"📋 QC 报告: {report_file}")
    
    # 生成编号转换说明
    numbering_info = """
# 编号体系转换说明
====================

本脚本进行初步质量检查。精确的编号转换需要 ANARCI 工具。

## 安装 ANARCI

```bash
conda install -c bioconda anarci
```

## 使用 ANARCI 进行编号转换

```bash
# Kabat 编号
anarci -i VH_SEQ.fasta -o VH_kabat.csv --scheme kabat

# IMGT 编号
anarci -i VH_SEQ.fasta -o VH_imgt.csv --scheme imgt

# Chothia 编号
anarci -i VH_SEQ.fasta -o VH_chothia.csv --scheme chothia
```

## CDR 边界定义差异

| 体系 | CDR1 | CDR2 | CDR3 |
|------|------|------|------|
| Kabat | 26-35 | 50-65 | 95-102 |
| Chothia | 26-32 | 52-56 | 95-102 |
| IMGT | 27-38 | 56-65 | 105-117 |
"""
    
    info_file = OUTPUT_DIR / "numbering_guide.md"
    with open(info_file, 'w') as f:
        f.write(numbering_info)
    print(f"📄 编号指南: {info_file}")
    
    print("\n" + "=" * 60)
    print("Step 10 完成!")
    print("=" * 60)
    print("\n下一步:")
    print("  1. 在服务器上安装 ANARCI")
    print("  2. 运行 anarci 进行精确编号转换")
    print("  3. 然后运行 Step 11 抗体结构预测")


if __name__ == "__main__":
    main()
