INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../shared/utils.cmake")
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/common.cmake")
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../shared/packages_nightly_performance.cmake")

set_cache_var(Trilinos_EXTRA_LINK_FLAGS "$ENV{CRAY_XPMEM_POST_LINK_OPTS} -lxpmem $ENV{PE_MPICH_GTL_DIR_amd_gfx90a} $ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}" STRING)

set_bool_cache_var(Trilinos_ENABLE_Fortran OFF)
set_cache_var(CMAKE_C_COMPILER "mpicc" STRING)
set_cache_var(CMAKE_CXX_COMPILER "mpicxx" STRING)

set_bool_cache_var(Tpetra_ENABLE_MMM_Timings ON)
set_bool_cache_var(Tpetra_ASSUME_GPU_AWARE_MPI ON)
set_bool_cache_var(Teuchos_TIMER_KOKKOS_FENCE ON)


set_bool_cache_var(Trilinos_ENABLE_OpenMP OFF)
set_bool_cache_var(Kokkos_ENABLE_OPENMP OFF)
set_cache_var(Phalanx_KOKKOS_DEVICE_TYPE "HIP" STRING)
set_bool_cache_var(Tpetra_INST_SERIAL ON)
set_bool_cache_var(Kokkos_ENABLE_HIP ON)
set_bool_cache_var(Tpetra_INST_HIP ON)
set_bool_cache_var(Sacado_ENABLE_HIERARCHICAL_DFAD ON)
set_bool_cache_var(TPL_ENABLE_ROCBLAS ON)
set_bool_cache_var(TPL_ENABLE_ROCSOLVER ON)
set_bool_cache_var(TPL_ENABLE_ROCSPARSE ON)


set_bool_cache_var(Kokkos_ARCH_ZEN3 ON)
set_bool_cache_var(Kokkos_ARCH_AMD_GFX90A ON)
