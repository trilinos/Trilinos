set_cache_var(CMAKE_CXX_FLAGS "-I$ENV{MPICH_DIR}/include" STRING)

set_bool_cache_var(TPL_ENABLE_MPI ON)
set_cache_var(MPI_EXEC "srun" STRING)
set_cache_var(MPI_EXEC_NUMPROCS_FLAG "-c7;--gpus-per-task=1;--gpu-bind=closest;-n" STRING)

set_bool_cache_var(TPL_ENABLE_BLAS ON)
set_cache_var(TPL_BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_amd.so" STRING)

set_bool_cache_var(TPL_ENABLE_LAPACK ON)
set_cache_var(TPL_LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_amd.so" STRING)

set_bool_cache_var(TPL_ENABLE_Netcdf ON)
set_cache_var(Netcdf_INCLUDE_DIRS "$ENV{NETCDF_DIR}/include;$ENV{PNETCDF_DIR}/include" STRING)
set_cache_var(Netcdf_LIBRARY_DIRS "$ENV{NETCDF_DIR}/lib;$ENV{PNETCDF_DIR}/lib" STRING)
set_bool_cache_var(TPL_Netcdf_PARALLEL ON)

set_bool_cache_var(TPL_ENABLE_Pnetcdf ON)
set_cache_var(Pnetcdf_INCLUDE_DIRS "$ENV{PNETCDF_DIR}/include" STRING)
set_cache_var(Pnetcdf_LIBRARY_DIRS "$ENV{PNETCDF_DIR}/lib" STRING)

set_bool_cache_var(TPL_ENABLE_HDF5 ON)
set_cache_var(HDF5_INCLUDE_DIRS "$ENV{HDF5_DIR}/include" STRING)
set_cache_var(HDF5_LIBRARY_DIRS "$ENV{HDF5_DIR}/lib" STRING)

set_bool_cache_var(TPL_ENABLE_gtest OFF)

set_bool_cache_var(BUILD_SHARED_LIBS OFF)
