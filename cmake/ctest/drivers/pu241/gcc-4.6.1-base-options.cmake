#
# Base options for all GCC 4.6.1 builds
#

# Define the core compilers
IF (NOT TRILINOS_TOOLSET_BASE)
  SET(TRILINOS_TOOLSET_BASE  /projects/vera/gcc-4.6.1/toolset)
ENDIF()
SET(GCC_BASE_DIR ${TRILINOS_TOOLSET_BASE}/gcc-4.6.1)
# Add rpath for compiler libraries and gomp for parts built with OpenMP
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS
  "-lgomp -Wl,-rpath,${GCC_BASE_DIR}/lib64"
  CACHE STRING "")
# This dicates downstream the intel fortran compiler to be used
# Include MKL and TBB; these should match version of Intel compilers being used
MESSAGE(STATUS "Selecting gfortran 4.6.1 compiler and Intel 12.0.4 TBB/MKL")
SET(IFORT_VERSION       "")
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-12.0.4-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-12.0.4-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with gfortran")
SET(LAPACK_INCLUDE_DIRS ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with gfortran")

# Build shared libs by default to save massive amounts of disk space
SET(BUILD_SHARED_LIBS ON CACHE BOOL
  "Set by default in gcc-4.6.1-base-options.cmake")

# To avoid problem with EpetraExt_inout_test failure in optimized code for hybrid builds
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL
  "Set by default in gcc-4.6.1-base-options.cmake")

# Turn off HDF5 in EpetraExt to avoid hdf5 conflicts
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL
  "Set by default in gcc-4.6.1-base-options.cmake")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
