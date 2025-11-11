#!/bin/bash
# 文件名: collect_loom_files.sh
# 功能：自动收集并重命名 Velocyto loom 文件

# 输入参数
input_root="03-cellranger"     # Cell Ranger 结果根目录
output_dir="05-velocyto-loom"  # 输出目录

# 创建输出目录
mkdir -p "$output_dir"

# 查找所有 SRX 样本目录
find "$input_root" -maxdepth 1 -type d -name "SRX*" | while read srx_dir; do
    # 提取 SRX ID
    srx_id=$(basename "$srx_dir")
    
    # 源文件路径
    src_file="${srx_dir}/10x/velocyto/10x.loom"
    
    # 目标路径
    dest_file="${output_dir}/${srx_id}.loom"
    
    # 执行复制
    if [[ -f "$src_file" ]]; then
        echo "正在处理: $srx_id"
        cp -v "$src_file" "$dest_file"
    else
        echo "⚠️ 文件不存在: $src_file"
    fi
done

echo "所有 loom 文件已保存至: $(realpath $output_dir)"
