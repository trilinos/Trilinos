if(TPL_ENABLE_MPI)
  set_cache_var(MPI_EXEC "srun" STRING)
  set_cache_var(MPI_EXEC_NUMPROCS_FLAG "-c32;--cpu-bind=cores;-G4;-n" STRING)
  set_cache_var(MPI_EXEC_POST_NUMPROCS_FLAG "${CMAKE_CURRENT_SOURCE_DIR}/../../environments/perlmutter/select_gpu_device.sh" STRING)

  set_cache_var(CMAKE_CXX_FLAGS "-I$ENV{MPICH_DIR}/include -I$CRAY_PE_LIBSCI_PREFIX/include" STRING)
  set_cache_var(CMAKE_EXE_LINKER_FLAGS "-L$ENV{MPICH_DIR}/lib -lmpi $ENV{PE_MPICH_GTL_DIR_nvidia80} $ENV{PE_MPICH_GTL_LIBS_nvidia80} -L$ENV{MPICH_DIR}/lib -lmpi_gnu_123 -L$ENV{MPICH_DIR}/lib/ -lmpichf90 -L$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_gnu -L$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_gnu_mpi" STRING)
endif()

set_cache_var(TPL_BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.so" STRING)

set_cache_var(TPL_LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.so" STRING)

set_cache_var(CMAKE_C_COMPILER "cc" STRING)
set_cache_var(CMAKE_CXX_COMPILER "$ENV{NVCC_WRAPPER}" STRING)
set_cache_var(CMAKE_Fortran_COMPILER "ftn" STRING)
