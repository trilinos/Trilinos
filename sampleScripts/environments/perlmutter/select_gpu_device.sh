#!/bin/bash
# select_gpu_device wrapper script
export CUDA_VISIBLE_DEVICES=$SLURM_LOCALID
exec $*
