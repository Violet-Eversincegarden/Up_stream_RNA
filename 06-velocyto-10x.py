#!/usr/bin/env python3
"""
Velocyto 自动化分析流程（含失败重试版）
功能：对失败样本使用 veloctyo run 命令重试，指定 BAM 和 barcode 路径
"""

import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def validate_paths(sample_dir, genes_gtf, repeat_gtf):
    """验证所有输入路径有效性"""
    checks = [
        (sample_dir, "Cell Ranger 输出目录"),
        (genes_gtf, "基因注释 GTF 文件"),
        (repeat_gtf, "重复序列 GTF 文件")
    ]
    
    missing = []
    for path, desc in checks:
        if not path.exists():
            missing.append(f"{desc}: {path}")
    
    if missing:
        print("❌ 路径验证失败:")
        print("\n".join(missing))
        return False
    return True

def run_velocyto_sample(args):
    """执行单个样本的 Velocyto 分析（使用 run10x）"""
    srx_dir, params = args
    print(f"\n🌀 开始处理 {srx_dir} (run10x)")
    
    # 配置路径参数
    cellranger_out = params['cellranger_root'] / srx_dir / "10x"
    output_dir = params['output_root'] / srx_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 验证关键路径
    if not validate_paths(cellranger_out, params['genes_gtf'], params['repeat_gtf']):
        return False
    
    # 构建 run10x 命令
    cmd = [
        "velocyto", "run10x",
        "--samtools-threads", str(params['samtools_threads']),
        "-v" if params['verbose'] else "",
        "-m", str(params['repeat_gtf']),
        str(cellranger_out),
        str(params['genes_gtf'])
    ]
    cmd = [c for c in cmd if c]  # 清理空参数
    
    # 运行命令
    log_file = output_dir / "velocyto_run10x.log"
    try:
        with open(log_file, 'w') as f:
            process = subprocess.run(
                cmd,
                check=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                cwd=output_dir,
                text=True
            )
        print(f"✅ {srx_dir} run10x 分析成功")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {srx_dir} run10x 分析失败 (code {e.returncode})，日志: {log_file}")
        return False

def retry_failed_sample(args):
    """使用 velocyto run 重试失败样本"""
    srx_dir, params = args
    print(f"\n🔄 重试失败样本 {srx_dir} (velocyto run)")
    
    # 配置路径
    cellranger_out = params['cellranger_root'] / srx_dir / "10x"
    output_dir = params['output_root'] / srx_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 准备文件路径
    bam_file = cellranger_out / "outs" / "possorted_genome_bam.bam"
    barcode_file = cellranger_out / "outs" / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    
    # 验证关键文件
    missing_files = []
    if not bam_file.exists():
        missing_files.append(f"BAM 文件: {bam_file}")
    if not barcode_file.exists():
        missing_files.append(f"Barcode 文件: {barcode_file}")
    
    if missing_files:
        print(f"❌ {srx_dir} 缺少必要文件:")
        print("\n".join(missing_files))
        return False
    
    # 构建 velocyto run 命令
    cmd = [
        "velocyto", "run",
        "--samtools-threads", str(params['samtools_threads']),
        "-v" if params['verbose'] else "",
        "-o", str(output_dir),
        "-b", str(barcode_file),
        "-m", str(params['repeat_gtf']),
        "-e", f"{srx_dir}",
        str(bam_file),
        str(params['genes_gtf'])
    ]
    cmd = [c for c in cmd if c]  # 清理空参数
    
    # 运行命令
    log_file = output_dir / "velocyto_retry.log"
    try:
        with open(log_file, 'w') as f:
            process = subprocess.run(
                cmd,
                check=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )
        print(f"✅ {srx_dir} 重试成功")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {srx_dir} 重试失败 (code {e.returncode})，日志: {log_file}")
        return False

def main():
    # 参数解析
    parser = argparse.ArgumentParser(description="Velocyto 并行分析工具（含失败重试）")
    parser.add_argument("-i", "--input", default="03-cellranger",
                      help="Cell Ranger 结果根目录（默认：03-cellranger）")
    parser.add_argument("-o", "--output", default="04-velocyto",
                      help="输出根目录（默认：04-velocyto）")
    parser.add_argument("-g", "--genes-gtf", required=True,
                      help="基因注释 GTF 文件路径（例：/path/to/genes.gtf）")
    parser.add_argument("-m", "--repeat-gtf", required=True,
                      help="重复序列 GTF 文件路径（例：/path/to/mm10_rmsk.gtf）")
    parser.add_argument("-j", "--max-workers", type=int, default=2,
                      help="最大并行任务数（默认：2）")
    parser.add_argument("-s", "--samtools-threads", type=int, default=2,
                      help="Samtools 线程数（默认：2）")
    parser.add_argument("-v", "--verbose", action="store_true",
                      help="启用详细输出")
    parser.add_argument("-r", "--retry-failed", action="store_true",
                      help="只重试之前失败的样本")
    args = parser.parse_args()

    # 初始化路径参数
    params = {
        'cellranger_root': Path(args.input).resolve(),
        'output_root': Path(args.output).resolve(),
        'genes_gtf': Path(args.genes_gtf).resolve(),
        'repeat_gtf': Path(args.repeat_gtf).resolve(),
        'samtools_threads': args.samtools_threads,
        'verbose': args.verbose
    }
    
    # 预验证全局路径
    global_checks = [
        (params['cellranger_root'], "Cell Ranger 输入目录"),
        (params['genes_gtf'], "基因注释文件"),
        (params['repeat_gtf'], "重复序列文件")
    ]
    for path, desc in global_checks:
        if not path.exists():
            print(f"❌ 全局路径错误: {desc} {path} 不存在")
            return

    # 获取样本列表
    srx_dirs = [d.name for d in params['cellranger_root'].iterdir() 
               if d.is_dir() and d.name.startswith("SRX")]
    print(f"🔍 发现 {len(srx_dirs)} 个样本")
    print(f"⚙️ 分析参数：")
    print(f" | Samtools 线程: {args.samtools_threads}")
    print(f" | 并行任务数: {args.max_workers}")
    print(f" | 详细模式: {'是' if args.verbose else '否'}")
    print(f" | 仅重试模式: {'是' if args.retry_failed else '否'}")

    # 模式选择
    if args.retry_failed:
        # 只重试之前失败的样本（输出目录存在但无 loom 文件）
        failed_samples = []
        for srx in srx_dirs:
            output_dir = params['output_root'] / srx
            loom_file = output_dir / f"{srx}.loom"
            if output_dir.exists() and not loom_file.exists():
                failed_samples.append(srx)
        
        print(f"\n🔧 准备重试 {len(failed_samples)} 个失败样本")
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            task_args = [(srx, params) for srx in failed_samples]
            results = list(executor.map(retry_failed_sample, task_args))
        
        retry_success = sum(results)
        print(f"\n{'='*40}")
        print(f"📊 重试完成统计:")
        print(f" ✅ 成功重试: {retry_success}")
        print(f" ❌ 重试失败: {len(failed_samples) - retry_success}")
    
    else:
        # 正常执行 + 自动重试
        # 首次运行所有样本
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            task_args = [(srx, params) for srx in srx_dirs]
            results = list(executor.map(run_velocyto_sample, task_args))

        # 识别失败样本
        failed_samples = [srx_dirs[i] for i, success in enumerate(results) if not success]
        success_count = sum(results)
        
        print(f"\n{'='*40}")
        print(f"📊 首次运行统计:")
        print(f" ✅ 成功: {success_count}")
        print(f" ❌ 失败: {len(failed_samples)}")
        
        # 如果有失败样本则重试
        if failed_samples:
            print(f"\n🔄 准备重试失败样本...")
            with ThreadPoolExecutor(max_workers=max(1, args.max_workers//2)) as executor:  # 减少重试线程数
                task_args = [(srx, params) for srx in failed_samples]
                retry_results = list(executor.map(retry_failed_sample, task_args))
            
            retry_success = sum(retry_results)
            success_count += retry_success
            print(f"\n{'='*40}")
            print(f"📊 最终统计:")
            print(f" ✅ 总计成功: {success_count}")
            print(f" ❌ 最终失败: {len(failed_samples) - retry_success}")

    print(f"📂 输出目录: {params['output_root']}")

if __name__ == "__main__":
    main()