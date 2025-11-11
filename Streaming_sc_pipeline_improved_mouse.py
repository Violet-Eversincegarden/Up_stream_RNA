#!/usr/bin/env python3
"""
改进版流式单细胞分析流程控制器
主要改进：
1. 添加依赖检查
2. 优化并发控制（分离下载和处理阶段）
3. 完善日志管理
4. 添加中间文件清理
5. 添加FastQC质控
6. 添加Cell Ranger指标检查
7. 支持断点续传
8. 动态资源分配
9. 更好的错误处理
"""

import os
import sys
import csv
import json
import argparse
import subprocess
import shutil
import time
import psutil
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# 默认配置
DEFAULT_CONFIG = {
    'threads': 16,
    'mem_gb': 64,
    'max_download_workers': 4,
    'max_processing_workers': 1,  # Cell Ranger很消耗资源，建议串行
    'transcriptome': "/public/home/tyj20231151/RNA-seq/ref.10x/refdata-gex-GRCm39-2024-A",
    'genes_gtf': "/public/home/tyj20231151/RNA-seq/ref.10x/refdata-gex-GRCm39-2024-A/genes/genes.gtf",
    'repeat_gtf': "/public/home/tyj20231151/RNA-seq/rmsk/mm10_rmsk.gtf",
    'prefetch_options': "--max-size u --progress",
    'enable_qc': True,
    'enable_cleanup': True,
    'checkpoint_enabled': True
}

# 全局变量
GLOBAL_CONFIG = DEFAULT_CONFIG.copy()
GLOBAL_STRATEGY = None

def setup_logging(base_dir):
    """设置日志目录"""
    log_dir = base_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建主日志文件
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    main_log = log_dir / f"pipeline_{timestamp}.log"
    
    return log_dir, main_log

def log_message(message, log_file=None, print_msg=True):
    """统一的日志记录函数"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] {message}"
    
    if print_msg:
        print(message)
    
    if log_file:
        with open(log_file, 'a') as f:
            f.write(log_line + '\n')

def check_dependencies():
    """检查必需工具是否安装"""
    required_tools = {
        'prefetch': 'SRA Toolkit',
        'fasterq-dump': 'SRA Toolkit',
        'pigz': 'Parallel gzip',
        'cellranger': 'Cell Ranger',
        'velocyto': 'Velocyto',
        'fastqc': 'FastQC (可选，用于质控)'
    }
    
    missing = []
    optional_missing = []
    
    print("\n🔍 检查依赖工具...")
    for tool, description in required_tools.items():
        if shutil.which(tool) is None:
            if tool == 'fastqc':
                optional_missing.append(f"{tool} ({description})")
            else:
                missing.append(f"{tool} ({description})")
            print(f"  ❌ {tool}: 未找到")
        else:
            print(f"  ✅ {tool}: 已安装")
    
    if missing:
        print(f"\n❌ 缺少必需工具: {', '.join(missing)}")
        print("请安装缺失的工具后重试")
        return False
    
    if optional_missing:
        print(f"\n⚠️ 缺少可选工具: {', '.join(optional_missing)}")
        print("将跳过相关功能")
    
    print("✅ 依赖检查通过\n")
    return True

def get_system_resources():
    """获取系统资源信息"""
    cpu_count = os.cpu_count() or 8
    mem_gb = psutil.virtual_memory().total / (1024**3)
    
    print(f"💻 系统资源:")
    print(f"  CPU核心数: {cpu_count}")
    print(f"  总内存: {mem_gb:.1f} GB")
    
    return cpu_count, mem_gb

def calculate_optimal_resources(cpu_count, mem_gb, max_workers):
    """计算最优的资源分配"""
    # 为每个任务分配资源，确保不超过系统容量
    per_task_threads = max(1, cpu_count // max_workers)
    per_task_mem = max(32, int(mem_gb * 0.8) // max_workers)  # 使用80%的内存
    
    print(f"📊 资源分配方案:")
    print(f"  并发任务数: {max_workers}")
    print(f"  每任务线程数: {per_task_threads}")
    print(f"  每任务内存: {per_task_mem} GB")
    
    return per_task_threads, per_task_mem

def determine_global_strategy(sra_groups, base_dir, log_file):
    """为整个GSE确定全局处理策略"""
    global GLOBAL_STRATEGY
    if GLOBAL_STRATEGY is not None:
        return GLOBAL_STRATEGY
    
    log_message("\n🔍 为整个GSE确定全局处理策略", log_file)
    
    # 找到第一个SRR用于策略测试
    test_srr = None
    for srx_id, srrs in sra_groups.items():
        if srrs:
            test_srr = srrs[0]
            break
    
    if not test_srr:
        log_message("⚠️ 没有可用的SRR进行策略测试", log_file)
        return None
    
    raw_dir = base_dir / "01-raw_fastq"
    raw_dir.mkdir(exist_ok=True)
    
    log_message(f"  使用SRR {test_srr} 进行策略测试...", log_file)
    
    # 尝试使用split-3策略
    cmd = f"fasterq-dump {test_srr} --outdir {raw_dir} --threads {GLOBAL_CONFIG['threads']} --split-3 --skip-technical --temp {raw_dir}"
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        generated_files = list(raw_dir.glob(f"{test_srr}*.fastq"))
        
        if len(generated_files) == 2:
            log_message(f"✅ {test_srr} 使用split-3产生2个文件，选择全局split-3策略", log_file)
            GLOBAL_STRATEGY = "split-3"
            
            # 清理测试文件
            for f in generated_files:
                f.unlink()
            return "split-3"
        else:
            for f in generated_files:
                f.unlink()
    except Exception as e:
        log_message(f"⚠️ {test_srr} split-3转换失败: {str(e)}", log_file)
    
    # 尝试使用split-files策略
    log_message(f"  尝试split-files策略...", log_file)
    try:
        cmd = f"fasterq-dump {test_srr} --outdir {raw_dir} --threads {GLOBAL_CONFIG['threads']} --split-files --include-technical --skip-technical --temp {raw_dir}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        generated_files = list(raw_dir.glob(f"{test_srr}*.fastq"))
        
        if len(generated_files) >= 2:
            log_message(f"✅ {test_srr} 使用split-files产生 {len(generated_files)} 个文件，选择全局split-files策略", log_file)
            GLOBAL_STRATEGY = "split-files"
            
            for f in generated_files:
                f.unlink()
            return "split-files"
        else:
            log_message(f"❌ {test_srr} split-files策略只产生 {len(generated_files)} 个文件", log_file)
    except Exception as e:
        log_message(f"❌ {test_srr} split-files策略转换失败: {str(e)}", log_file)
    
    return None

class SRXProcessor:
    """处理单个SRX的完整流程"""
    
    def __init__(self, srx_id, srrs, base_dir, config, log_dir):
        self.srx_id = srx_id
        self.srrs = srrs
        self.file_types = set()
        self.base_dir = Path(base_dir).resolve()
        self.config = config
        self.log_dir = log_dir
        
        # 创建SRX专用日志文件
        self.log_file = self.log_dir / f"{srx_id}.log"
        
        # 定义目录
        self.raw_dir = self.base_dir / "01-raw_fastq"
        self.merged_dir = self.base_dir / "02-merged_fastq" / srx_id
        self.qc_dir = self.base_dir / "02-fastqc" / srx_id
        self.cellranger_dir = self.base_dir / "03-cellranger" / srx_id
        self.matrix_dir = self.base_dir / "04-filtered_matrices"
        self.velocyto_dir = self.base_dir / "04-velocyto" / srx_id
        self.loom_dir = self.base_dir / "05-velocyto-loom"
        
        # 创建必要的目录
        for dir_path in [self.raw_dir, self.merged_dir, self.qc_dir, 
                         self.cellranger_dir, self.matrix_dir, 
                         self.velocyto_dir, self.loom_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # 状态追踪（用于断点续传）
        self.checkpoint_file = self.log_dir / f"{srx_id}_checkpoint.json"
        self.status = self.load_checkpoint()
    
    def log(self, message, level="INFO"):
        """记录日志"""
        log_message(f"[{self.srx_id}] {message}", self.log_file)
    
    def load_checkpoint(self):
        """加载检查点"""
        if self.config.get('checkpoint_enabled') and self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    status = json.load(f)
                self.log(f"📋 从检查点恢复: {sum(status.values())}/{len(status)} 步骤已完成")
                return status
            except:
                pass
        
        return {
            'download': False,
            'convert': False,
            'merge': False,
            'qc': False,
            'cellranger': False,
            'check_metrics': False,
            'extract_matrix': False,
            'velocyto': False,
            'collect_loom': False,
            'cleanup': False
        }
    
    def save_checkpoint(self):
        """保存检查点"""
        if self.config.get('checkpoint_enabled'):
            with open(self.checkpoint_file, 'w') as f:
                json.dump(self.status, f, indent=2)
    
    def _run_command(self, cmd, log_name, shell=True):
        """统一的命令执行接口"""
        cmd_log = self.log_dir / self.srx_id / log_name
        cmd_log.parent.mkdir(parents=True, exist_ok=True)
        
        self.log(f"执行命令: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
        
        with open(cmd_log, 'w') as f:
            result = subprocess.run(cmd, shell=shell, stdout=f, stderr=subprocess.STDOUT)
        
        if result.returncode != 0:
            self.log(f"命令失败，详见日志: {cmd_log}")
        
        return result.returncode == 0
    
    def download_srrs(self):
        """下载SRR文件"""
        if self.status['download']:
            self.log("⏭️ 跳过已完成的步骤: 下载")
            return True
        
        self.log(f"⬇️ 开始下载 {len(self.srrs)} 个SRR文件")
        
        for srr in self.srrs:
            cmd = f"prefetch {srr} {self.config['prefetch_options']}"
            if not self._run_command(cmd, f"prefetch_{srr}.log"):
                self.log(f"❌ {srr} 下载失败")
                return False
            self.log(f"✅ {srr} 下载完成")
        
        self.status['download'] = True
        self.save_checkpoint()
        return True
    
    def convert_fastq(self):
        """转换SRR为FASTQ"""
        if self.status['convert']:
            self.log("⏭️ 跳过已完成的步骤: FASTQ转换")
            return True
        
        global GLOBAL_STRATEGY
        self.log(f"🔄 转换FASTQ文件（策略: {GLOBAL_STRATEGY}）")
        
        for srr in self.srrs:
            try:
                if GLOBAL_STRATEGY == "split-3":
                    cmd = f"fasterq-dump {srr} --outdir {self.raw_dir} --threads {self.config['threads']} --split-3 --skip-technical --temp {self.raw_dir}"
                    if not self._run_command(cmd, f"fasterq_{srr}.log"):
                        return False
                    
                    generated_files = list(self.raw_dir.glob(f"{srr}*.fastq"))
                    if len(generated_files) != 2:
                        self.log(f"❌ {srr} 产生了 {len(generated_files)} 个文件，期望2个")
                        return False
                    
                    # 重命名
                    file1 = self.raw_dir / f"{srr}_1.fastq"
                    file2 = self.raw_dir / f"{srr}_2.fastq"
                    
                    if not (file1.exists() and file2.exists()):
                        self.log(f"❌ {srr} 缺少必需文件")
                        return False
                    
                    file1.rename(self.raw_dir / f"{srr}_R1.fastq")
                    file2.rename(self.raw_dir / f"{srr}_R2.fastq")
                    
                    # 压缩
                    subprocess.run(["pigz", "-f", "-p", str(self.config['threads']), 
                                  str(self.raw_dir / f"{srr}_R1.fastq")])
                    subprocess.run(["pigz", "-f", "-p", str(self.config['threads']), 
                                  str(self.raw_dir / f"{srr}_R2.fastq")])
                    
                    self.file_types.update(["R1", "R2"])
                
                elif GLOBAL_STRATEGY == "split-files":
                    cmd = f"fasterq-dump {srr} --outdir {self.raw_dir} --threads {self.config['threads']} --split-files --include-technical --skip-technical --temp {self.raw_dir}"
                    if not self._run_command(cmd, f"fasterq_{srr}.log"):
                        return False
                    
                    generated_files = list(self.raw_dir.glob(f"{srr}*.fastq"))
                    if not generated_files:
                        self.log(f"❌ {srr} 未生成任何文件")
                        return False
                    
                    mapping = self.determine_file_mapping(srr)
                    if not mapping:
                        return False
                    
                    for suffix, name in mapping.items():
                        src = self.raw_dir / f"{srr}_{suffix}.fastq"
                        dest = self.raw_dir / f"{srr}_{name}.fastq"
                        
                        if src.exists():
                            src.rename(dest)
                            if os.path.getsize(dest) > 0:
                                subprocess.run(["pigz", "-f", "-p", str(self.config['threads']), str(dest)])
                            else:
                                self.log(f"⚠️ 文件 {dest.name} 为空，跳过压缩")
                
                self.log(f"✅ {srr} 转换完成")
            
            except Exception as e:
                self.log(f"❌ {srr} 转换失败: {str(e)}")
                return False
        
        self.status['convert'] = True
        self.save_checkpoint()
        return True
    
    def determine_file_mapping(self, srr):
        """根据实际文件数量动态确定文件类型映射"""
        file_patterns = [f"{srr}_1.fastq", f"{srr}_2.fastq", 
                        f"{srr}_3.fastq", f"{srr}_4.fastq"]
        existing_files = [f for f in file_patterns if (self.raw_dir / f).exists()]
        num_files = len(existing_files)
        
        if num_files == 3:
            mapping = {"1": "I1", "2": "R1", "3": "R2"}
        elif num_files == 4:
            seq_lengths = {}
            for num in ["1", "2", "3", "4"]:
                file_path = self.raw_dir / f"{srr}_{num}.fastq"
                if file_path.exists():
                    try:
                        with open(file_path, 'r') as f:
                            for i, line in enumerate(f):
                                if i % 4 == 1:
                                    seq_lengths[num] = len(line.strip())
                                    break
                    except Exception:
                        continue
            
            if len(seq_lengths) < 2:
                self.log(f"❌ {srr} 无法确定有效文件")
                return None
            
            sorted_nums = sorted(seq_lengths.keys(), 
                               key=lambda x: seq_lengths.get(x, 0), reverse=True)[:2]
            mapping = {sorted_nums[0]: "R1", sorted_nums[1]: "R2"}
            
            # 删除不使用的文件
            for num in ["1", "2", "3", "4"]:
                if num not in sorted_nums:
                    file_path = self.raw_dir / f"{srr}_{num}.fastq"
                    if file_path.exists():
                        file_path.unlink()
                        self.log(f"🗑️ 删除非主要读段: {file_path.name}")
        
        elif num_files == 2:
            mapping = {"1": "R1", "2": "R2"}
        elif num_files == 1:
            mapping = {"1": "R1"}
        else:
            self.log(f"❌ 不支持的文件数量: {num_files}")
            return None
        
        self.file_types.update(mapping.values())
        self.log(f"  文件映射: {mapping}")
        return mapping
    
    def merge_fastq(self):
        """合并同一SRX的FASTQ文件"""
        if self.status['merge']:
            self.log("⏭️ 跳过已完成的步骤: FASTQ合并")
            return True
        
        self.log(f"🧩 合并FASTQ文件")
        
        file_types = self.file_types if self.file_types else ["R1", "R2"]
        
        for file_type in file_types:
            output_file = self.merged_dir / f"{self.srx_id}_S1_L001_{file_type}_001.fastq.gz"
            with open(output_file, "wb") as out_f:
                files_found = False
                
                for srr in self.srrs:
                    src_file = self.raw_dir / f"{srr}_{file_type}.fastq.gz"
                    
                    if src_file.exists():
                        files_found = True
                        with open(src_file, "rb") as in_f:
                            shutil.copyfileobj(in_f, out_f)
                        self.log(f"  ➕ 添加 {src_file.name}")
                
                if not files_found:
                    self.log(f"⚠️ 未找到 {file_type} 类型文件")
                    try:
                        output_file.unlink()
                    except:
                        pass
                else:
                    self.log(f"  ✅ 合并 {file_type} 完成")
        
        self.status['merge'] = True
        self.save_checkpoint()
        return True
    
    def run_fastqc(self):
        """运行FastQC质量控制"""
        if self.status['qc']:
            self.log("⏭️ 跳过已完成的步骤: 质控")
            return True
        
        if not self.config.get('enable_qc', True):
            self.log("⏭️ 质控已禁用，跳过")
            self.status['qc'] = True
            return True
        
        if shutil.which('fastqc') is None:
            self.log("⚠️ FastQC未安装，跳过质控")
            self.status['qc'] = True
            return True
        
        self.log("🔍 运行FastQC质量控制")
        
        fastq_files = list(self.merged_dir.glob("*.fastq.gz"))
        if not fastq_files:
            self.log("⚠️ 未找到FASTQ文件，跳过质控")
            self.status['qc'] = True
            return True
        
        for fastq in fastq_files:
            cmd = f"fastqc {fastq} -o {self.qc_dir} -t {self.config['threads']}"
            if not self._run_command(cmd, f"fastqc_{fastq.stem}.log"):
                self.log(f"⚠️ FastQC失败: {fastq.name}")
        
        self.log(f"✅ 质控完成，结果保存在: {self.qc_dir}")
        self.status['qc'] = True
        self.save_checkpoint()
        return True
    
    def run_cellranger(self):
        """运行Cell Ranger分析"""
        if self.status['cellranger']:
            self.log("⏭️ 跳过已完成的步骤: Cell Ranger")
            return True
        
        self.log(f"🔬 运行Cell Ranger分析")
        
        if not any(os.scandir(self.merged_dir)):
            self.log(f"❌ FASTQ目录为空: {self.merged_dir}")
            return False
        
        cmd = [
            "cellranger", "count",
            "--id", "10x",
            "--transcriptome", self.config['transcriptome'],
            "--fastqs", str(self.merged_dir),
            "--sample", self.srx_id,
            "--create-bam", "true",
            "--localcores", str(self.config['threads']),
            "--localmem", str(self.config['mem_gb'])
        ]
        
        if not self._run_command(cmd, "cellranger.log", shell=False):
            self.log("❌ Cell Ranger分析失败")
            return False
        
        self.log("✅ Cell Ranger分析完成")
        self.status['cellranger'] = True
        self.save_checkpoint()
        return True
    
    def check_cellranger_metrics(self):
        """检查Cell Ranger输出指标"""
        if self.status['check_metrics']:
            self.log("⏭️ 跳过已完成的步骤: 指标检查")
            return True
        
        self.log("📊 检查Cell Ranger输出指标")
        
        metrics_file = self.cellranger_dir / "10x" / "outs" / "metrics_summary.csv"
        if not metrics_file.exists():
            self.log("⚠️ 未找到metrics_summary.csv")
            self.status['check_metrics'] = True
            return True
        
        try:
            with open(metrics_file, 'r') as f:
                reader = csv.DictReader(f)
                metrics = next(reader)
            
            # 提取关键指标
            key_metrics = {
                'Estimated Number of Cells': metrics.get('Estimated Number of Cells', 'N/A'),
                'Mean Reads per Cell': metrics.get('Mean Reads per Cell', 'N/A'),
                'Median Genes per Cell': metrics.get('Median Genes per Cell', 'N/A'),
                'Total Genes Detected': metrics.get('Total Genes Detected', 'N/A'),
                'Reads Mapped Confidently to Transcriptome': metrics.get('Reads Mapped Confidently to Transcriptome', 'N/A')
            }
            
            self.log("  关键指标:")
            for key, value in key_metrics.items():
                self.log(f"    {key}: {value}")
            
            # 简单的质量检查
            try:
                n_cells = int(metrics.get('Estimated Number of Cells', '0').replace(',', ''))
                if n_cells < 100:
                    self.log(f"⚠️ 警告: 细胞数量较少 ({n_cells})")
            except:
                pass
            
        except Exception as e:
            self.log(f"⚠️ 解析指标文件失败: {str(e)}")
        
        self.status['check_metrics'] = True
        self.save_checkpoint()
        return True
    
    def extract_matrix(self):
        """提取feature矩阵"""
        if self.status['extract_matrix']:
            self.log("⏭️ 跳过已完成的步骤: 矩阵提取")
            return True
        
        self.log(f"📊 提取特征矩阵")
        
        possible_paths = [
            self.cellranger_dir / "10x" / "outs" / "filtered_feature_bc_matrix",
            self.cellranger_dir / "10x" / "outs" / "filtered_gene_bc_matrices"
        ]
        
        src_dir = None
        for path in possible_paths:
            if path.exists():
                src_dir = path
                break
        
        if not src_dir:
            self.log(f"❌ 未找到矩阵目录")
            return False
        
        dest_dir = self.matrix_dir / f"{self.srx_id}_matrix"
        
        try:
            if dest_dir.exists():
                shutil.rmtree(dest_dir)
            shutil.copytree(src_dir, dest_dir)
            self.log(f"✅ 矩阵提取完成: {dest_dir}")
            
            self.status['extract_matrix'] = True
            self.save_checkpoint()
            return True
        except Exception as e:
            self.log(f"❌ 矩阵提取失败: {str(e)}")
            return False
    
    def run_velocyto(self):
        """运行Velocyto分析"""
        if self.status['velocyto']:
            self.log("⏭️ 跳过已完成的步骤: Velocyto")
            return True
        
        self.log(f"⚡ 运行Velocyto分析")
        
        bam_file = self.cellranger_dir / "10x" / "outs" / "possorted_genome_bam.bam"
        if not bam_file.exists():
            self.log(f"❌ 缺少BAM文件: {bam_file}")
            return False
        
        cmd = [
            "velocyto", "run10x",
            "--samtools-threads", str(self.config['threads']),
            "-m", self.config['repeat_gtf'],
            str(self.cellranger_dir / "10x"),
            self.config['genes_gtf']
        ]
        
        # Velocyto需要在cellranger目录下运行
        if not self._run_command(cmd, "velocyto.log", shell=False):
            self.log("❌ Velocyto分析失败")
            return False
        
        self.log(f"✅ Velocyto分析完成")
        self.status['velocyto'] = True
        self.save_checkpoint()
        return True
    
    def collect_loom(self):
        """收集loom文件"""
        if self.status['collect_loom']:
            self.log("⏭️ 跳过已完成的步骤: loom收集")
            return True
        
        self.log(f"📦 收集loom文件")
        
        src_file = self.cellranger_dir / "10x" / "velocyto" / "10x.loom"
        dest_file = self.loom_dir / f"{self.srx_id}.loom"
        
        if not src_file.exists():
            self.log(f"❌ 未找到loom文件: {src_file}")
            return False
        
        try:
            if dest_file.exists():
                dest_file.unlink()
            shutil.copy(src_file, dest_file)
            self.log(f"✅ 已复制 {src_file.name} -> {dest_file}")
            
            self.status['collect_loom'] = True
            self.save_checkpoint()
            return True
        except Exception as e:
            self.log(f"❌ loom文件收集失败: {str(e)}")
            return False
    
    def cleanup(self):
        """清理中间文件"""
        if self.status['cleanup']:
            self.log("⏭️ 跳过已完成的步骤: 清理")
            return True
        
        if not self.config.get('enable_cleanup', True):
            self.log("⏭️ 清理已禁用，跳过")
            self.status['cleanup'] = True
            return True
        
        self.log("🧹 清理中间文件")
        
        cleaned_size = 0
        
        # 清理未压缩的FASTQ
        for f in self.raw_dir.glob(f"*{self.srx_id}*.fastq"):
            size = f.stat().st_size
            f.unlink()
            cleaned_size += size
            self.log(f"  🗑️ 删除: {f.name}")
        
        # 清理SRA文件
        sra_dir = Path.home() / "ncbi/public/sra"
        if sra_dir.exists():
            for srr in self.srrs:
                sra_file = sra_dir / f"{srr}.sra"
                if sra_file.exists():
                    size = sra_file.stat().st_size
                    sra_file.unlink()
                    cleaned_size += size
                    self.log(f"  🗑️ 删除: {sra_file.name}")
        
        # 清理原始FASTQ压缩文件（如果不再需要）
        for srr in self.srrs:
            for file_type in self.file_types:
                raw_gz = self.raw_dir / f"{srr}_{file_type}.fastq.gz"
                if raw_gz.exists():
                    size = raw_gz.stat().st_size
                    raw_gz.unlink()
                    cleaned_size += size
                    self.log(f"  🗑️ 删除: {raw_gz.name}")
        
        cleaned_gb = cleaned_size / (1024**3)
        self.log(f"✅ 清理完成，释放空间: {cleaned_gb:.2f} GB")
        
        self.status['cleanup'] = True
        self.save_checkpoint()
        return True
    
    def process_full_workflow(self):
        """执行完整的处理流程"""
        try:
            start_time = time.time()
            self.log(f"\n{'='*60}")
            self.log(f"🚀 开始处理 (包含 {len(self.srrs)} 个SRR)")
            self.log(f"  工作目录: {self.base_dir}")
            self.log(f"{'-'*60}")
            
            # 定义处理步骤
            steps = [
                ('download', '下载', self.download_srrs),
                ('convert', 'FASTQ转换', self.convert_fastq),
                ('merge', 'FASTQ合并', self.merge_fastq),
                ('qc', '质量控制', self.run_fastqc),
                ('cellranger', 'Cell Ranger分析', self.run_cellranger),
                ('check_metrics', '指标检查', self.check_cellranger_metrics),
                ('extract_matrix', '矩阵提取', self.extract_matrix),
                ('velocyto', 'Velocyto分析', self.run_velocyto),
                ('collect_loom', 'loom收集', self.collect_loom),
                ('cleanup', '清理', self.cleanup)
            ]
            
            for step_id, step_name, step_func in steps:
                try:
                    if not step_func():
                        self.log(f"❌ {step_name}失败")
                        return False
                except Exception as e:
                    self.log(f"❌ {step_name}异常: {str(e)}")
                    return False
            
            duration = time.time() - start_time
            self.log(f"{'-'*60}")
            self.log(f"✅ 处理完成! 耗时: {duration/60:.1f}分钟")
            self.log(f"{'='*60}\n")
            
            return True
        
        except Exception as e:
            self.log(f"❌ 处理失败: {str(e)}")
            return False

def load_metadata(metadata_file):
    """加载元数据并分组SRX"""
    srx_groups = defaultdict(list)
    try:
        with open(metadata_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                srx = row["Experiment"]
                srr = row["Run"]
                srx_groups[srx].append(srr)
        return srx_groups
    except Exception as e:
        print(f"❌ 加载元数据失败: {str(e)}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="改进版流式单细胞分析流程控制器",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s -m SraRunTable.csv
  %(prog)s -m SraRunTable.csv -t 16 -g 64
  %(prog)s -m SraRunTable.csv --max-processing-workers 2
  %(prog)s -m SraRunTable.csv --disable-qc --disable-cleanup
        """
    )
    
    parser.add_argument("-m", "--metadata", default="SraRunTable.csv",
                       help="元数据文件 (CSV格式)")
    parser.add_argument("-t", "--threads", type=int,
                       help=f"每个任务的CPU线程数 (默认: {DEFAULT_CONFIG['threads']})")
    parser.add_argument("-g", "--mem-gb", type=int,
                       help=f"每个任务的内存(GB) (默认: {DEFAULT_CONFIG['mem_gb']})")
    parser.add_argument("-d", "--max-download-workers", type=int,
                       help=f"下载并发数 (默认: {DEFAULT_CONFIG['max_download_workers']})")
    parser.add_argument("-p", "--max-processing-workers", type=int,
                       help=f"处理并发数 (默认: {DEFAULT_CONFIG['max_processing_workers']})")
    parser.add_argument("--disable-qc", action="store_true",
                       help="禁用FastQC质量控制")
    parser.add_argument("--disable-cleanup", action="store_true",
                       help="禁用中间文件清理")
    parser.add_argument("--disable-checkpoint", action="store_true",
                       help="禁用断点续传")
    parser.add_argument("--auto-resources", action="store_true",
                       help="自动计算最优资源分配")
    
    args = parser.parse_args()
    
    # 检查依赖
    if not check_dependencies():
        sys.exit(1)
    
    # 更新全局配置
    global GLOBAL_CONFIG
    if args.threads:
        GLOBAL_CONFIG['threads'] = args.threads
    if args.mem_gb:
        GLOBAL_CONFIG['mem_gb'] = args.mem_gb
    if args.max_download_workers:
        GLOBAL_CONFIG['max_download_workers'] = args.max_download_workers
    if args.max_processing_workers:
        GLOBAL_CONFIG['max_processing_workers'] = args.max_processing_workers
    if args.disable_qc:
        GLOBAL_CONFIG['enable_qc'] = False
    if args.disable_cleanup:
        GLOBAL_CONFIG['enable_cleanup'] = False
    if args.disable_checkpoint:
        GLOBAL_CONFIG['checkpoint_enabled'] = False
    
    # 自动资源分配
    if args.auto_resources:
        cpu_count, mem_gb = get_system_resources()
        per_task_threads, per_task_mem = calculate_optimal_resources(
            cpu_count, mem_gb, GLOBAL_CONFIG['max_processing_workers']
        )
        GLOBAL_CONFIG['threads'] = per_task_threads
        GLOBAL_CONFIG['mem_gb'] = per_task_mem
    
    # 获取工作目录
    base_dir = Path.cwd().resolve()
    
    # 设置日志
    log_dir, main_log = setup_logging(base_dir)
    
    print("="*60)
    print(f"🔬 启动改进版流式单细胞分析流程")
    print(f"⚙️ 配置参数:")
    print(f"  线程数: {GLOBAL_CONFIG['threads']}")
    print(f"  内存(GB): {GLOBAL_CONFIG['mem_gb']}")
    print(f"  下载并发数: {GLOBAL_CONFIG['max_download_workers']}")
    print(f"  处理并发数: {GLOBAL_CONFIG['max_processing_workers']}")
    print(f"  质量控制: {'启用' if GLOBAL_CONFIG['enable_qc'] else '禁用'}")
    print(f"  文件清理: {'启用' if GLOBAL_CONFIG['enable_cleanup'] else '禁用'}")
    print(f"  断点续传: {'启用' if GLOBAL_CONFIG['checkpoint_enabled'] else '禁用'}")
    print(f"📁 工作目录: {base_dir}")
    print(f"📁 元数据文件: {args.metadata}")
    print(f"📁 日志目录: {log_dir}")
    print("="*60)
    
    # 加载元数据
    srx_groups = load_metadata(args.metadata)
    if srx_groups is None:
        sys.exit(1)
    
    total_srr = sum(len(v) for v in srx_groups.values())
    log_message(f"🔍 发现 {len(srx_groups)} 个SRX组，包含 {total_srr} 个SRR", main_log)
    
    # 确定全局策略
    global_strategy = determine_global_strategy(srx_groups, base_dir, main_log)
    if not global_strategy:
        log_message("❌ 无法确定全局策略，退出流程", main_log)
        sys.exit(1)
    
    # 阶段1: 批量下载
    log_message("\n" + "="*60, main_log)
    log_message("📥 阶段1: 批量下载SRR文件", main_log)
    log_message("="*60, main_log)
    
    processors = []
    for srx_id, srrs in srx_groups.items():
        processor = SRXProcessor(srx_id, srrs, base_dir, GLOBAL_CONFIG, log_dir)
        processors.append(processor)
    
    download_success = 0
    with ThreadPoolExecutor(max_workers=GLOBAL_CONFIG['max_download_workers']) as executor:
        futures = {executor.submit(p.download_srrs): p for p in processors}
        
        for future in as_completed(futures):
            processor = futures[future]
            try:
                if future.result():
                    download_success += 1
                    log_message(f"✅ {processor.srx_id} 下载完成", main_log)
                else:
                    log_message(f"❌ {processor.srx_id} 下载失败", main_log)
            except Exception as e:
                log_message(f"❌ {processor.srx_id} 下载异常: {str(e)}", main_log)
    
    log_message(f"\n📊 下载阶段完成: {download_success}/{len(processors)} 成功", main_log)
    
    if download_success == 0:
        log_message("❌ 所有下载失败，退出流程", main_log)
        sys.exit(1)
    
    # 阶段2: 分析处理
    log_message("\n" + "="*60, main_log)
    log_message("🔬 阶段2: 分析处理", main_log)
    log_message("="*60, main_log)
    
    success_count = 0
    with ThreadPoolExecutor(max_workers=GLOBAL_CONFIG['max_processing_workers']) as executor:
        futures = {}
        for processor in processors:
            # 跳过下载失败的
            if not processor.status['download']:
                continue
            future = executor.submit(processor.process_full_workflow)
            futures[future] = processor.srx_id
        
        for future in as_completed(futures):
            srx_id = futures[future]
            try:
                if future.result():
                    success_count += 1
                    log_message(f"✅ {srx_id} 流程完成", main_log)
                else:
                    log_message(f"❌ {srx_id} 流程失败", main_log)
            except Exception as e:
                log_message(f"❌ {srx_id} 发生异常: {str(e)}", main_log)
    
    # 最终报告
    print("\n" + "="*60)
    print(f"📊 流程汇总统计")
    print(f"✅ 成功: {success_count}/{len(srx_groups)}")
    print(f"❌ 失败: {len(srx_groups) - success_count}")
    print(f"\n📂 输出目录结构:")
    print(f"  {base_dir}/01-raw_fastq/       # 原始FASTQ文件")
    print(f"  {base_dir}/02-merged_fastq/    # 合并后的FASTQ")
    print(f"  {base_dir}/02-fastqc/          # FastQC质控报告")
    print(f"  {base_dir}/03-cellranger/      # Cell Ranger结果")
    print(f"  {base_dir}/04-filtered_matrices/ # 特征矩阵")
    print(f"  {base_dir}/04-velocyto/        # Velocyto结果")
    print(f"  {base_dir}/05-velocyto-loom/   # 收集的loom文件")
    print(f"  {base_dir}/logs/               # 日志文件")
    print("="*60)
    
    if success_count == len(srx_groups):
        print("🎉 所有SRX处理完成!")
        log_message("🎉 所有SRX处理完成!", main_log)
    else:
        print("⚠️ 部分SRX处理失败，请检查日志")
        log_message("⚠️ 部分SRX处理失败", main_log)
    
    print(f"\n📝 详细日志: {main_log}")

if __name__ == "__main__":
    main()