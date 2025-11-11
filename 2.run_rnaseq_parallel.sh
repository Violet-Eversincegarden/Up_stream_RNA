#!/bin/bash

# ==============================================================================
# RNA-seq 上游分析流程脚本 - 并发版本
#
# 新增功能:
#   - 支持多样本并发处理
#   - 可控制最大并发数
#   - 更好的日志管理
# ==============================================================================

set -euo pipefail

# ========================= 用户配置 =========================

BASE_DIR="/home/tyj/project/HQJ/01-raw-data" 
STAR_INDEX="/home/tyj/RNA-seq/ref.genecode.mouse/STAR"
ANNOTATION_GTF="/home/tyj/RNA-seq/ref.genecode.mouse/gencode.vM35.annotation.gtf" 
OUTPUT_DIR="/home/tyj/project/HQJ/analysis_results"
THREADS=32
STRANDEDNESS=0

# --- 新增并发控制参数 ---
# 最大同时处理的样本数（建议根据内存和CPU核心数设置）
MAX_PARALLEL_SAMPLES=4
# 每个样本使用的线程数（总线程数 = MAX_PARALLEL_SAMPLES * THREADS_PER_SAMPLE）
THREADS_PER_SAMPLE=$((THREADS / MAX_PARALLEL_SAMPLES))

# 确保每个样本至少有4个线程
if [ $THREADS_PER_SAMPLE -lt 4 ]; then
    THREADS_PER_SAMPLE=4
    MAX_PARALLEL_SAMPLES=$((THREADS / 4))
    echo "[WARNING] 调整并发参数: MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES}, THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE}"
fi

# ====================== 脚本主体 ======================

echo "======================================================"
echo "RNA-seq 并发分析流程开始"
echo "最大并发样本数: ${MAX_PARALLEL_SAMPLES}"
echo "每样本线程数: ${THREADS_PER_SAMPLE}"
echo "总计算线程数: $((MAX_PARALLEL_SAMPLES * THREADS_PER_SAMPLE))"
echo "======================================================"

# 前置检查
for tool in fastqc fastp STAR featureCounts multiqc parallel; do
    if ! command -v "$tool" &> /dev/null; then
        echo "[ERROR] $tool 未安装或不在PATH中"
        exit 1
    fi
done

# 输入验证
for path in "${BASE_DIR}" "${STAR_INDEX}"; do
    [ ! -d "$path" ] && { echo "[ERROR] 目录不存在: $path"; exit 1; }
done
[ ! -f "${ANNOTATION_GTF}" ] && { echo "[ERROR] 文件不存在: ${ANNOTATION_GTF}"; exit 1; }

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"/{01_fastqc_raw,02_fastp_trimmed,03_fastqc_trimmed,04_star_alignment,05_featurecounts,06_multiqc,logs}

# 定义单个样本处理函数
process_sample() {
    local SAMPLE_DIR="$1"
    local SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
    local LOG_FILE="${OUTPUT_DIR}/logs/${SAMPLE_NAME}.log"
    
    # 重定向该样本的所有输出到独立日志文件
    exec 1> >(tee -a "${LOG_FILE}")
    exec 2> >(tee -a "${LOG_FILE}" >&2)
    
    echo "[START] $(date): 开始处理样本 ${SAMPLE_NAME}"
    
    # 查找FASTQ文件
    mapfile -t READ1_FILES < <(find "${SAMPLE_DIR}" -maxdepth 1 -name "*_R1*.fastq.gz" -o -name "*_1.fastq.gz")
    mapfile -t READ2_FILES < <(find "${SAMPLE_DIR}" -maxdepth 1 -name "*_R2*.fastq.gz" -o -name "*_2.fastq.gz")
    
    if [ ${#READ1_FILES[@]} -ne 1 ] || [ ${#READ2_FILES[@]} -ne 1 ]; then
        echo "[ERROR] ${SAMPLE_NAME}: 未找到唯一的成对FASTQ文件"
        return 1
    fi
    
    local READ1="${READ1_FILES[0]}"
    local READ2="${READ2_FILES[0]}"
    
    # 验证文件配对
    local base1=$(basename "${READ1}" | sed 's/_R[12].*//g' | sed 's/_[12].*//g')
    local base2=$(basename "${READ2}" | sed 's/_R[12].*//g' | sed 's/_[12].*//g')
    if [ "${base1}" != "${base2}" ]; then
        echo "[ERROR] ${SAMPLE_NAME}: 文件名不匹配"
        return 1
    fi
    
    echo "[INFO] ${SAMPLE_NAME}: 找到配对文件"
    echo "  R1: ${READ1}"
    echo "  R2: ${READ2}"
    
    # 1. FastQC (原始数据)
    echo "[STEP 1/4] ${SAMPLE_NAME}: FastQC (原始数据)"
    fastqc -o "${OUTPUT_DIR}/01_fastqc_raw" -t "${THREADS_PER_SAMPLE}" "${READ1}" "${READ2}" || return 1
    
    # 2. fastp 清洗
    echo "[STEP 2/4] ${SAMPLE_NAME}: fastp 数据清洗"
    local TRIMMED_R1="${OUTPUT_DIR}/02_fastp_trimmed/${SAMPLE_NAME}_R1.trimmed.fastq.gz"
    local TRIMMED_R2="${OUTPUT_DIR}/02_fastp_trimmed/${SAMPLE_NAME}_R2.trimmed.fastq.gz"
    
    fastp -i "${READ1}" -I "${READ2}" \
          -o "${TRIMMED_R1}" -O "${TRIMMED_R2}" \
          --html "${OUTPUT_DIR}/02_fastp_trimmed/${SAMPLE_NAME}.fastp.html" \
          --json "${OUTPUT_DIR}/02_fastp_trimmed/${SAMPLE_NAME}.fastp.json" \
          -w "${THREADS_PER_SAMPLE}" || return 1
    
    # 3. FastQC (清洗后数据)
    echo "[INFO] ${SAMPLE_NAME}: FastQC (清洗后数据)"
    fastqc -o "${OUTPUT_DIR}/03_fastqc_trimmed" -t "${THREADS_PER_SAMPLE}" "${TRIMMED_R1}" "${TRIMMED_R2}" || true
    
    # 4. STAR 比对
    echo "[STEP 3/4] ${SAMPLE_NAME}: STAR 比对"
    STAR --genomeDir "${STAR_INDEX}" \
         --runThreadN "${THREADS_PER_SAMPLE}" \
         --readFilesIn "${TRIMMED_R1}" "${TRIMMED_R2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --quantMode GeneCounts || return 1
    
    # 5. featureCounts 定量
    echo "[STEP 4/4] ${SAMPLE_NAME}: featureCounts 定量"
    local BAM_FILE="${OUTPUT_DIR}/04_star_alignment/${SAMPLE_NAME}_Aligned.sortedByCoord.out.bam"
    
    featureCounts -T "${THREADS_PER_SAMPLE}" \
                  -p -s "${STRANDEDNESS}" \
                  -a "${ANNOTATION_GTF}" \
                  -o "${OUTPUT_DIR}/05_featurecounts/${SAMPLE_NAME}_counts.txt" \
                  "${BAM_FILE}" || return 1
    
    echo "[SUCCESS] $(date): 样本 ${SAMPLE_NAME} 处理完成"
    return 0
}

# 导出函数和变量供 parallel 使用
export -f process_sample
export OUTPUT_DIR STAR_INDEX ANNOTATION_GTF THREADS_PER_SAMPLE STRANDEDNESS

# 获取所有样本目录
mapfile -t SAMPLE_DIRS < <(find "${BASE_DIR}" -maxdepth 1 -type d -not -path "${BASE_DIR}")

if [ ${#SAMPLE_DIRS[@]} -eq 0 ]; then
    echo "[ERROR] 在 ${BASE_DIR} 中未找到任何样本目录"
    exit 1
fi

echo "[INFO] 找到 ${#SAMPLE_DIRS[@]} 个样本目录，开始并发处理..."

# 使用 GNU Parallel 并发处理
printf "%s\n" "${SAMPLE_DIRS[@]}" | \
parallel -j "${MAX_PARALLEL_SAMPLES}" \
         --bar \
         --joblog "${OUTPUT_DIR}/logs/parallel_joblog.txt" \
         --results "${OUTPUT_DIR}/logs/parallel_results" \
         process_sample {}

# 检查处理结果
TOTAL_SAMPLES=${#SAMPLE_DIRS[@]}
SUCCESS_COUNT=$(grep -c "SUCCESS" "${OUTPUT_DIR}"/logs/*.log 2>/dev/null || echo 0)
FAILED_COUNT=$((TOTAL_SAMPLES - SUCCESS_COUNT))

echo "======================================================"
echo "并发处理完成！"
echo "  - 总样本数: ${TOTAL_SAMPLES}"
echo "  - 成功处理: ${SUCCESS_COUNT}"
echo "  - 失败数量: ${FAILED_COUNT}"
echo "======================================================"


