#!/bin/bash
#SBATCH --job-name=HMM  # 作业名
#SBATCH --partition=gpu3090         # gpu3090队列
#SBATCH --qos=1gpu                  # gpu3090Qos
#SBATCH --nodes=1                   # 节点数量
#SBATCH --ntasks-per-node=1         # 每节点进程数
#SBATCH --cpus-per-task=4           # 1:4 的 GPU:CPU 配比 
#SBATCH --gres=gpu:1                # 2 块 GPU
#SBATCH --output=logs/%j.out        # 标准输出
#SBATCH --error=logs/%j.err         # 错误输出
#SBATCH --mail-user=yueqing.xing22@student.xjtlu.edu.cn
#SBATCH --mail-type=ALL

module load hmmer/3.3.2-gcc-8.5.0-7vk4bd2
module load mafft/7.505-gcc-8.5.0-omqojgg

PROJECT_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/scripts"
UNKNOWN_FASTA="/gpfs/work/bio/yueqingxing22/igem/HMM/data/unknown_sequences/proteins_converted.fasta"
HMM_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/hmm_profiles/pfam_families/hmm"
OUT_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/annotation_results/pfam_families"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "$PROJECT_DIR"
python functional_annotator.py annotate \
    --unknown_fasta "$UNKNOWN_FASTA" \
    --hmm_dir "$HMM_DIR" \
    --out_dir "$OUT_DIR" \
    --threads "$THREADS" \
    --iE 1e-5 --bits 25 --qcov 0.30
