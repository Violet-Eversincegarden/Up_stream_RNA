#!/usr/bin/env python3
"""
10x 数据分析脚本 (最终版)
功能：执行 Cell Ranger 分析并正确组织输出目录
"""

import csv
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def load_srx_mapping(metadata_file):
    """加载元数据"""
    srx_map = defaultdict(list)
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            srx = row["Experiment"]
            srr = row["Run"]
            srx_map[srx].append(srr)
    return srx_map

def run_cellranger(srx, input_root, output_root, transcriptome, cores, mem):
    """运行单个样本分析"""
    # 配置路径
    sample_input_dir = input_root / srx
    sample_output_dir = output_root / srx  # 注意，只建到SRX，不加10x
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    
    # 构建命令
    cmd = [
        "cellranger", "count",
        "--id", "10x",  # 这个ID对应的是 cellranger自己新建的子目录
        "--transcriptome", str(transcriptome.resolve()),
        "--fastqs", str(sample_input_dir.resolve()),  # 指向fastq目录
        "--sample", srx,
        "--create-bam", "true",
        "--localcores", str(cores),
        "--localmem", str(mem)
    ]
    
    # 运行命令，把工作目录设到 SRX 目录（不是10x目录！）
    log_file = sample_output_dir / "cellranger.log"
    try:
        with open(log_file, 'w') as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT, cwd=sample_output_dir, text=True)
        print(f"✅ {srx} 分析成功")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {srx} 分析失败，查看日志: {log_file}")
        print(f"错误详情: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="10x Cell Ranger 自动化工具")
    parser.add_argument("-m", "--metadata", default="SraRunTable.csv",
                        help="元数据文件")
    parser.add_argument("-i", "--input-dir", default="02-merged_fastq",
                        help="预处理后的输入目录")
    parser.add_argument("-o", "--output-root", default="03-cellranger",
                        help="分析结果输出根目录")
    parser.add_argument("-r", "--transcriptome", required=True,
                        help="参考转录组路径")
    parser.add_argument("-j", "--max-workers", type=int, default=2,
                        help="并行任务数")
    parser.add_argument("-c", "--per-task-cores", type=int, default=32,
                        help="每个任务使用的CPU核心数")
    parser.add_argument("-M", "--per-task-mem", type=int, default=128,
                        help="每个任务使用的内存(GB)")
    args = parser.parse_args()
    
    # 校验输入路径
    input_dir = Path(args.input_dir)
    transcriptome = Path(args.transcriptome)
    if not input_dir.exists():
        print(f"❌ 输入目录不存在: {input_dir}")
        return
    if not transcriptome.exists():
        print(f"❌ 参考转录组不存在: {transcriptome}")
        return

    # 加载元数据
    srx_map = load_srx_mapping(args.metadata)
    print(f"🔍 发现样本数: {len(srx_map)}")
    print(f"⚙️ 并行配置: {args.max_workers} 并发 | 每任务 {args.per_task_cores} 核 {args.per_task_mem}GB 内存")
    
    # 并行执行
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = []
        for srx in srx_map:
            future = executor.submit(
                run_cellranger,
                srx,
                input_dir,
                Path(args.output_root),
                transcriptome,
                args.per_task_cores,
                args.per_task_mem
            )
            futures.append(future)
        
        # 结果统计
        success_count = sum(f.result() for f in futures)
        print(f"\n🎉 分析完成: 成功 {success_count}/{len(srx_map)} 个样本")
        print(f"✅ 结果保存目录: {args.output_root}")

if __name__ == "__main__":
    main()
