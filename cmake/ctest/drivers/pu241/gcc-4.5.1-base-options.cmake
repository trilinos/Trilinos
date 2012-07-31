#
# Base options for all hybrid GCC 4.5.1 C/C++ and Intel 11.1.064 Fortran builds
#

# Define the core compilers
SET(TRILINOS_TOOLSET_BASE  /opt/gcc-4.5.1/trilinos-toolset)
# Add rpath for compiler libraries
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS "-Wl,-rpath,${TRILINOS_TOOLSET_BASE}/lib64" CACHE STRING "")
# This dicates downstream the intel fortran compiler to be used
# Include MKL and TBB; these should match version of Intel compilers being used
MESSAGE(STATUS "Selecting ifort 11.1.064 compiler and libraries")
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-11.1.064-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-11.1.064-options.cmake)
SET(IFORT_VERSION       "11.1.064")
SET(BLAS_INCLUDE_DIRS   ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

# To avoid problem with EpetraExt_inout_test failure in optimized code for hybrid build
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
