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
FAMILIES_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/data/families/pfam_families"
HMM_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/hmm_profiles/pfam_families/hmm"
ALN_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM/hmm_profiles/pfam_families/alignments"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "$PROJECT_DIR"
python functional_annotator.py build\
  --families_dir "$FAMILIES_DIR" \
  --aln_dir "$ALN_DIR" \
  --hmm_dir "$HMM_DIR" \
  --threads "$THREADS"