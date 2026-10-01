#!/usr/bin/env python3
"""
Phase 2 Pipeline 入口脚本
===========================
一键运行所有步骤

用法:
    python run_pipeline.py --from-step 10 --to-step 16
    python run_pipeline.py --step 11
    python run_pipeline.py --mode blind --to-step 13
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime


# ============================================================
# 配置
# ============================================================

SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent


# ============================================================
# 步骤定义
# ============================================================

STEPS = {
    10: {
        'name': '编号体系统一',
        'script': 'step10_numbering.py',
        'description': '质量检查 + ANARCI 编号转换'
    },
    11: {
        'name': '抗体结构预测',
        'script': 'step11_antibody_structure.py',
        'description': 'ABodyBuilder2 预测 VH/VL Fv 结构'
    },
    12: {
        'name': '复合物预测',
        'script': 'step12_complex_prediction.py',
        'description': 'AlphaFold3/Boltz-1 预测 IL-24+Ab 复合物'
    },
    13: {
        'name': '接触位点定位',
        'script': 'step13_contact_mapping.py',
        'description': 'd < 4.5 Å 距离判定 epitope/paratope'
    },
    14: {
        'name': '复合物质量评估',
        'script': 'step14_quality_assessment.py',
        'description': 'ipTM/PAE/ensemble 一致性'
    },
    15: {
        'name': '亲和力评估',
        'script': 'step15_affinity_estimation.py',
        'description': 'PRODIGY/Rosetta 能量估算'
    },
    16: {
        'name': '竞争性结合分析',
        'script': 'step16_competition_analysis.py',
        'description': '中和机制判定'
    }
}


# ============================================================
# 主函数
# ============================================================

def run_step(step_num: int, mode: str = 'blind', **kwargs):
    """运行单个步骤"""
    step_info = STEPS.get(step_num)
    if not step_info:
        print(f"❌ 未知步骤: {step_num}")
        return False
    
    script_path = SCRIPT_DIR / step_info['script']
    if not script_path.exists():
        print(f"❌ 脚本不存在: {script_path}")
        return False
    
    print("\n" + "=" * 60)
    print(f"Step {step_num}: {step_info['name']}")
    print("=" * 60)
    print(f"📄 {step_info['description']}")
    
    # 构建命令
    cmd = [sys.executable, str(script_path)]
    
    if step_num >= 12:
        cmd.extend(['--mode', f'12{mode[0]}'])  # 12A_blind 或 12B_guided
    
    # 运行
    result = os.system(' '.join(cmd))
    
    if result == 0:
        print(f"✅ Step {step_num} 完成")
        return True
    else:
        print(f"❌ Step {step_num} 失败 (exit code: {result})")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Phase 2 Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run_pipeline.py                          # 运行所有步骤
  python run_pipeline.py --from-step 11           # 从 Step 11 开始
  python run_pipeline.py --to-step 13             # 运行到 Step 13
  python run_pipeline.py --step 11                # 只运行 Step 11
  python run_pipeline.py --mode guided            # 使用 Guided 模式 (默认 blind)
        """
    )
    
    parser.add_argument('--from-step', type=int, default=10,
                       help='起始步骤 (默认: 10)')
    parser.add_argument('--to-step', type=int, default=16,
                       help='结束步骤 (默认: 16)')
    parser.add_argument('--step', type=int, default=None,
                       help='只运行指定步骤')
    parser.add_argument('--mode', type=str, default='blind',
                       choices=['blind', 'guided'],
                       help='复合物预测模式 (默认: blind)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Phase 2 Pipeline")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    # 确定要运行的步骤
    if args.step:
        steps_to_run = [args.step]
    else:
        steps_to_run = list(range(args.from_step, args.to_step + 1))
    
    print(f"📋 步骤: {steps_to_run}")
    print(f"🎯 模式: {args.mode}")
    
    # 运行每个步骤
    failed_steps = []
    for step_num in steps_to_run:
        success = run_step(step_num, mode=args.mode)
        if not success:
            failed_steps.append(step_num)
    
    # 总结
    print("\n" + "=" * 60)
    print("运行完成")
    print("=" * 60)
    print(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if failed_steps:
        print(f"❌ 失败步骤: {failed_steps}")
    else:
        print("✅ 所有步骤成功完成")
    
    # 输出目录
    print(f"\n📁 输出目录: {PROJECT_DIR / 'results'}")
    
    # 列出输出文件
    results_dir = PROJECT_DIR / "results"
    if results_dir.exists():
        print("\n📄 生成的文件:")
        for step_dir in sorted(results_dir.iterdir()):
            if step_dir.is_dir():
                files = list(step_dir.glob("*"))
                print(f"   {step_dir.name}/: {len(files)} 个文件")


if __name__ == "__main__":
    main()
