#!/usr/bin/env python3
"""
Velocyto 自动化分析流程（参数修复版）
功能：严格遵循 CLI 参数顺序，优化错误处理
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
    """执行单个样本的 Velocyto 分析"""
    srx_dir, params = args
    print(f"\n🌀 开始处理 {srx_dir}")
    
    # 配置路径参数
    cellranger_out = params['cellranger_root'] / srx_dir /"10x"
    output_dir = params['output_root'] / srx_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 验证关键路径
    if not validate_paths(cellranger_out, params['genes_gtf'], params['repeat_gtf']):
        return False
    
    # 构建命令（严格遵循参数顺序）
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
    log_file = output_dir / "velocyto.log"
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
        print(f"✅ {srx_dir} 分析成功")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {srx_dir} 分析失败 (code {e.returncode})，日志: {log_file}")
        return False

def main():
    # 参数解析
    parser = argparse.ArgumentParser(description="Velocyto 并行分析工具（修复版）")
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

    # 并行执行
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        task_args = [(srx, params) for srx in srx_dirs]
        results = list(executor.map(run_velocyto_sample, task_args))

    # 统计结果
    success_count = sum(results)
    print(f"\n{'='*40}")
    print(f"📊 任务完成统计:")
    print(f" ✅ 成功: {success_count}")
    print(f" ❌ 失败: {len(srx_dirs) - success_count}")
    print(f"📂 输出目录: {params['output_root']}")

if __name__ == "__main__":
    main()