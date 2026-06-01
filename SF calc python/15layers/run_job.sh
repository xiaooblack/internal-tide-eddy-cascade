#!/bin/bash


# 设置并行线程数（关键！）
export NUMBA_NUM_THREADS=64
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64

# ============ 打印作业信息 ============
echo "=========================================="
echo "作业ID: $SLURM_JOB_ID"
echo "作业名: $SLURM_JOB_NAME"
echo "节点: $SLURM_NODELIST"
echo "CPU核心数: $SLURM_CPUS_PER_TASK"
echo "开始时间: $(date)"
echo "工作目录: $(pwd)"
echo "=========================================="

# ============ 运行程序 ============
python -u main_15layers.py 0 69

# ============ 结束 ============
echo "=========================================="
echo "结束时间: $(date)"
echo "=========================================="
