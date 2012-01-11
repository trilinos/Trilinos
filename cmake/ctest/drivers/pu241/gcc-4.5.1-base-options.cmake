#
# Base options for all hybrid GCC 4.5.1 C/C++ and Intel Fortran builds
#

# Define the core compilers
SET(TRILINOS_TOOLSET_BASE  /opt/gcc-4.5.1/trilinos-toolset)
# Add rpath for compiler libraries
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS "-Wl,-rpath,${TRILINOS_TOOLSET_BASE}/lib64" CACHE STRING "")
# This dicates downstream the intel fortran compiler to be used
# Include MKL and TBB; these should match version of Intel compilers being used
IF ("${HYBRIDBUILD_INTEL_VERSION}" STREQUAL "")
  MESSAGE(STATUS "Selecting gfortran 4.5.1 compiler and Intel 12.0.4 TBB/MKL")
  SET(IFORT_VERSION       "")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-12.0.4-options.cmake)
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-12.0.4-options.cmake)
  SET(BLAS_INCLUDE_DIRS   ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with gfortran")
  SET(LAPACK_INCLUDE_DIRS ${MKL_GCC451_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with gfortran")
ELSE()
  IF ("${HYBRIDBUILD_INTEL_VERSION}" STREQUAL "12.0.4")
    MESSAGE(STATUS "Selecting ifort 12.0.4 compiler and libraries")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-12.0.4-options.cmake)
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-12.0.4-options.cmake)
  ELSEIF("${HYBRIDBUILD_INTEL_VERSION}" STREQUAL "11.1.064")
    MESSAGE(STATUS "Selecting ifort 11.1.064 compiler and libraries")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-11.1.064-options.cmake)
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-11.1.064-options.cmake)
  ENDIF()
  SET(IFORT_VERSION       ${HYBRIDBUILD_INTEL_VERSION})
  SET(BLAS_INCLUDE_DIRS   ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with Intel fortran")
  SET(LAPACK_INCLUDE_DIRS ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with Intel fortran")
ENDIF()

# To avoid problem with EpetraExt_inout_test failure in optimized code for hybrid build
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL "")

# this compiler supports BinUtils
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")

# Point to basic CASL-related TPLs related to the GCC C/C++ 4.5.1 compiler
SET(PVMLibraries_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/pvm3/lib/LINUX64 CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/hdf5-1.8.5-patch1/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/hdf5-1.8.5-patch1/include CACHE FILEPATH "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-vri-tpls.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
