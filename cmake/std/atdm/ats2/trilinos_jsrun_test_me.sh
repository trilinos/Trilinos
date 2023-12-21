#!/bin/bash

echo "testing with gpu-aware unset!, should not see -gpu, look for -disable_gpu_hooks"
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun -M -disable_gdr -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun -M -disable_gdr -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun -M -disable_gdr -p 1 hostname

echo "testing with gpu-aware turned off!, should not see -gpu, look for -disable_gpu_hooks"
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun -M -disable_gdr -p 1 hostname

echo "testing with gpu-aware turned ON!, should see -gpu if NP>1 and  -disable_gpu_hooks when NP=1"
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun -M -disable_gdr -p 1 hostname

echo "testing with gpu-aware unset!, should not see -gpu, look for -disable_gpu_hooks"
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun  -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun  -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI="" ./trilinos_jsrun  -p 1 hostname

echo "testing with gpu-aware turned off!, should not see -gpu, look for -disable_gpu_hooks"
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun  -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun  -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=0 ./trilinos_jsrun  -p 1 hostname

echo "testing with gpu-aware turned ON!, should see -gpu if NP>1 and  -disable_gpu_hooks when NP=1"
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun  -p 2 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun  -p 1 hostname
TPETRA_ASSUME_GPU_AWARE_MPI=1 ./trilinos_jsrun  -p 1 hostname
