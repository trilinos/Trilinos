INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../shared/utils.cmake")

set_bool_cache_var(BUILD_SHARED_LIBS OFF)
set_bool_cache_var(Trilinos_ENABLE_Fortran OFF)

set_bool_cache_var(Kokkos_ARCH_ZEN3 ON)
set_bool_cache_var(Kokkos_ARCH_AMPERE80 ON)
set_bool_cache_var(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC OFF)

set_bool_cache_var(Tpetra_ENABLE_MMM_Timings ON)
set_bool_cache_var(Tpetra_ASSUME_GPU_AWARE_MPI ON)
set_bool_cache_var(Teuchos_TIMER_KOKKOS_FENCE ON)

set_bool_cache_var(Sacado_ENABLE_HIERARCHICAL_DFAD ON)

set_bool_cache_var(TPL_ENABLE_CUDA ON)
set_bool_cache_var(TPL_ENABLE_CUBLAS ON)
set_bool_cache_var(TPL_ENABLE_CUSPARSE ON)
set_bool_cache_var(TPL_ENABLE_CUSOLVER ON)
set_bool_cache_var(TPL_ENABLE_MPI ON)
set_bool_cache_var(TPL_ENABLE_BLAS ON)
set_bool_cache_var(TPL_ENABLE_LAPACK ON)
set_bool_cache_var(TPL_ENABLE_HDF5 ON)
set_bool_cache_var(TPL_ENABLE_Netcdf ON)
set_bool_cache_var(TPL_ENABLE_gtest OFF)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/common.cmake")
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../shared/packages_nightly_performance.cmake")
