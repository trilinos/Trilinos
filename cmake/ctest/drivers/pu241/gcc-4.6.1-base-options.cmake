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
# The ANC/VIPRE/BOA code does not work with gfortran 4.6.1 (or any GCC version of Fortran)
SET(VERA_ENABLE_CASLRAVE OFF CACHE BOOL "")
SET(VERA_ENABLE_CASLBOA OFF CACHE BOOL "")

# To avoid problem with EpetraExt_inout_test failure in optimized code for hybrid builds
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL "")

# This compiler supports BinUtils
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")

# Point to basic CASL-related TPLs related to the GCC C/C++ 4.6.1 compiler
# SET(PVMLibraries_LIBRARY_DIRS /opt/gcc-4.6.1/tpls/pvm3/lib/LINUX64 CACHE FILEPATH "")
# SET(PVMHeaders_INCLUDE_DIRS /opt/gcc-4.6.1/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/gcc-4.6.1/tpls/hdf5-1.8.9/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/gcc-4.6.1/tpls/hdf5-1.8.9/include CACHE FILEPATH "")
SET(Netcdf_INCLUDE_DIRS /opt/gcc-4.6.1/tpls/netcdf-4.2/include CACHE FILEPATH "")
SET(Netcdf_LIBRARY_DIRS /opt/gcc-4.6.1/tpls/netcdf-4.2/lib     CACHE FILEPATH "")
SET(Zlib_INCLUDE_DIRS   /opt/gcc-4.6.1/tpls/zlib-1.2.7/include CACHE FILEPATH "")
SET(Zlib_LIBRARY_DIRS   /opt/gcc-4.6.1/tpls/zlib-1.2.7/lib     CACHE FILEPATH "")
SET(QT_REQUIRED_VERSION 4.7.1                           CACHE STRING   "")
SET(QT_QMAKE_EXECUTABLE /opt/gcc-4.6.1/tpls/qt-4.7.1/bin/qmake CACHE FILEPATH "")
SET(MOOSE_PETSC_INCLUDE_DIRS  /opt/gcc-4.6.1/tpls/petsc-3.1-p8/include  CACHE FILEPATH "")
SET(MOOSE_PETSC_LIBRARY_DIRS  /opt/gcc-4.6.1/tpls/petsc-3.1-p8/lib      CACHE FILEPATH "")
SET(MOOSE_HYPRE_INCLUDE_DIRS  /opt/gcc-4.6.1/tpls/hypre-2.8.0b/include  CACHE FILEPATH "")
SET(MOOSE_HYPRE_LIBRARY_DIRS  /opt/gcc-4.6.1/tpls/hypre-2.8.0b/lib      CACHE FILEPATH "")
# SET(QT_LIBRARY_DIR /opt/gcc-4.6.1/tpls/qt-4.7.1/lib     CACHE FILEPATH "")
# SET(QT_INCLUDE_DIR /opt/gcc-4.6.1/tpls/qt-4.7.1/include CACHE FILEPATH "")

# Used only when TriKota/Dakota are enabled with CMake
SET(TriKota_ENABLE_DakotaCMake ON CACHE BOOL "")
SET(DAKOTA_ENABLE_TESTS OFF CACHE BOOL "")
SET(PECOS_ENABLE_TESTS OFF CACHE BOOL "")
SET(LHS_ENABLE_TESTS OFF CACHE BOOL "")
SET(HOPSPACK_ENABLE_TESTS OFF CACHE BOOL "")
SET(OPTPP_ENABLE_TESTS OFF CACHE BOOL "")
SET(HAVE_ACRO OFF CACHE BOOL "")
SET(HAVE_AMPL OFF CACHE BOOL "")
SET(HAVE_X_GRAPHICS OFF CACHE BOOL "")
SET(HAVE_HOPSPACK OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-vri-tpls.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
