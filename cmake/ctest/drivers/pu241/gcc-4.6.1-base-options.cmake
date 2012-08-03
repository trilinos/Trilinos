#
# Base options for all GCC 4.6.1 builds
#

# Define the core compilers
SET(TRILINOS_TOOLSET_BASE  /opt/gcc-4.6.1/trilinos-toolset)
# Add rpath for compiler libraries and gomp for parts built with OpenMP
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS
  "-lgomp -Wl,-rpath,${TRILINOS_TOOLSET_BASE}/lib64"
  CACHE STRING "")
# This dicates downstream the intel fortran compiler to be used
# Include MKL and TBB; these should match version of Intel compilers being used
MESSAGE(STATUS "Selecting gfortran 4.6.1 compiler and Intel 12.0.4 TBB/MKL")
SET(IFORT_VERSION       "")
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-12.0.4-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-12.0.4-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with gfortran")
SET(LAPACK_INCLUDE_DIRS ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with gfortran")

# To avoid problem with EpetraExt_inout_test failure in optimized code for hybrid builds
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL "")

# Turn off HDF5 in EpetraExt to avoid hdf5 conflicts
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
