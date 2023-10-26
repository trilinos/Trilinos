module load sparc-dev/cuda-11.2.152_gcc-8.3.1_spmpi-rolling
# sparc-dev sets NVCC_WRAPPER_TMPDIR to /tmp/${USER}
# This is what we want (avoids possible filename conflicts with other users),
# but just need to make sure the directory exists.
mkdir -p /tmp/${USER}
export ATDM_CONFIG_SYSTEM_NAME=ats2
export CUDA_LAUNCH_BLOCKING=0
export TPETRA_ASSUME_GPU_AWARE_MPI=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank

# Temporary!
# When #12447 is merged, this can be safely deleted.
# Until then, it's needed for trilinos_jsrun to pass the correct
# "-M -gpu" flag to mpiexec.
export TPETRA_ASSUME_CUDA_AWARE_MPI=1
