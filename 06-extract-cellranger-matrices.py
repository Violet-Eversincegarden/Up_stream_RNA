#!/usr/bin/env python3
"""
提取所有 Cell Ranger 的 filtered_feature_bc_matrix 到统一目录 (SRX 专用版)
"""

import csv
import shutil
from pathlib import Path
import argparse
import os

def extract_matrices(metadata_file, cellranger_output, output_root):
    """根据 SRX ID 提取所有矩阵数据"""
    # 1. 从元数据文件中读取所有 SRX ID
    srx_list = []
    with open(metadata_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            srx = row["Experiment"]
            if srx not in srx_list:
                srx_list.append(srx)
    
    print(f"🔍 发现 {len(srx_list)} 个 SRX 样本")
    
    # 2. 准备输出目录
    output_dir = Path(output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    extracted_count = 0
    failed_count = 0
    
    # 3. 对于每个 SRX 样本
    for srx in srx_list:
        # 3.1 构建矩阵路径 (使用 SRX 作为目录名)
        matrix_path = Path(cellranger_output) / srx / "10x" / "outs" / "filtered_feature_bc_matrix"
        
        if not matrix_path.exists():
            print(f"❌ 找不到 {srx} 的矩阵数据: {matrix_path}")
            failed_count += 1
            continue
        
        # 3.2 目标路径
        dest_path = output_dir / f"{srx}_matrix"
        
        try:
            # 删除可能存在的旧目录
            if dest_path.exists():
                shutil.rmtree(dest_path)
            
            # 复制整个矩阵目录
            shutil.copytree(matrix_path, dest_path)
            print(f"✅ 提取成功: {srx} -> {dest_path}")
            
            # 验证基本文件存在
            required_files = {"barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"}
            existing_files = set(f.name for f in dest_path.iterdir())
            
            if not required_files.issubset(existing_files):
                print(f"⚠️ 警告: {srx} 的矩阵数据可能不完整")
                failed_count += 1
                # 删除不完整数据
                shutil.rmtree(dest_path)
            else:
                extracted_count += 1
        except Exception as e:
            print(f"❌ 复制失败 {srx}: {str(e)}")
            failed_count += 1
    
    print("\n" + "="*50)
    print(f"成功提取: {extracted_count} 个样本")
    print(f"失败: {failed_count} 个样本")
    print(f"输出目录: {output_dir.resolve()}")
    
    # 生成样本列表文件
    sample_list = output_dir / "sample_list.txt"
    with open(sample_list, "w") as f:
        for srx in srx_list:
            sample_dir = output_dir / f"{srx}_matrix"
            if sample_dir.exists():
                f.write(f"{srx}\t{sample_dir}\n")
    
    print(f"样本列表已保存: {sample_list}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="提取 Cell Ranger 分析结果矩阵")
    parser.add_argument("-m", "--metadata", default="SraRunTable.csv", 
                        help="元数据文件 (必须包含 Experiment 列)")
    parser.add_argument("-i", "--cellranger-output", default="03-cellranger", 
                        help="Cell Ranger 输出根目录")
    parser.add_argument("-o", "--output-root", default="04-filtered_matrices", 
                        help="矩阵输出目录")
    
    args = parser.parse_args()
    extract_matrices(
        metadata_file=args.metadata,
        cellranger_output=args.cellranger_output,
        output_root=args.output_root
    )