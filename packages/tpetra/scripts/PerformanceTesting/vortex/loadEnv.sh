module load sparc-dev/cuda-11.2.152_gcc-8.3.1_spmpi-rolling
export NVCC_WRAPPER_TMPDIR=/tmp
export ATDM_CONFIG_SYSTEM_NAME=ats2
export CUDA_LAUNCH_BLOCKING=0
export TPETRA_ASSUME_GPU_AWARE_MPI=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
