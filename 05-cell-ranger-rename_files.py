#!/usr/bin/env python3
"""
10x 数据预处理脚本 (重命名与移动文件)
功能：创建标准目录结构并重命名文件（不可逆操作）
"""

import csv
import argparse
import shutil
from pathlib import Path
from collections import defaultdict

def load_srx_mapping(metadata_file):
    """加载元数据，构建 SRX 到 SRR 列表的映射"""
    srx_map = defaultdict(list)
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            srx = row["Experiment"]
            srr = row["Run"]
            srx_map[srx].append(srr)
    return srx_map

def prepare_sample_dirs(srx, merged_dir):
    """为每个 SRX 创建目录并重命名文件"""
    # 创建目标目录
    srx_dir = merged_dir / srx
    srx_dir.mkdir(exist_ok=True)
    
    # 定义文件映射规则
    file_map = {
        f"{srx}_1.fastq.gz": f"{srx}_S1_L001_R1_001.fastq.gz",
        f"{srx}_2.fastq.gz": f"{srx}_S1_L001_R2_001.fastq.gz"
    }
    
    # 移动并重命名文件
    for src_name, dest_name in file_map.items():
        src_path = merged_dir / src_name
        dest_path = srx_dir / dest_name
        
        if not src_path.exists():
            print(f"❌ 文件缺失: {src_path}")
            return False
        
        if dest_path.exists():
            print(f"⚠️ 目标文件已存在: {dest_path}，跳过")
            continue
            
        shutil.move(str(src_path), str(dest_path))
        print(f"📂 移动文件: {src_path} → {dest_path}")
    
    return True

def main():
    # 参数解析
    parser = argparse.ArgumentParser(description="10x 数据预处理工具")
    parser.add_argument("-m", "--metadata", default="SraRunTable.csv",
                      help="元数据文件（默认：SraRunTable.csv）")
    parser.add_argument("-i", "--input-dir", default="02-merged_fastq",
                      help="输入文件目录（默认：02-merged_fastq）")
    args = parser.parse_args()
    
    # 校验输入目录
    merged_dir = Path(args.input_dir)
    if not merged_dir.exists():
        print(f"❌ 输入目录不存在: {merged_dir}")
        return
    
    # 加载元数据
    srx_map = load_srx_mapping(args.metadata)
    print(f"🔍 发现 {len(srx_map)} 个 SRX 样本")
    
    # 处理每个样本
    for srx in srx_map:
        print(f"\n{'='*40}")
        print(f"🔄 处理样本: {srx}")
        success = prepare_sample_dirs(srx, merged_dir)
        print(f"✅ 完成" if success else "❌ 失败")

if __name__ == "__main__":
    main()