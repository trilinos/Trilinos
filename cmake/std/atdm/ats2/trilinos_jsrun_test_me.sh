#!/bin/bash

echo "testing with cuda-aware turned off!, should not see -gpu, look for -disable_gpu_hooks"
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 2 hostname 1>/dev/null
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 1 hostname 1>/dev/null
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 1 hostname 1>/dev/null

echo "testing with cuda-aware turned ON!, should see -gpu if NP>1 and  -disable_gpu_hooks when NP=1"
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 2 hostname 1>/dev/null
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 1 hostname 1>/dev/null
ECHO_CMD=1 TPETRA_ASSUME_CUDA_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 1 hostname 1>/dev/null
