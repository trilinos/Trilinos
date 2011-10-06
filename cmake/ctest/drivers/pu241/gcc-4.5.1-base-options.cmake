
#
# Base options for all hybrid GCC 4.5.1 C/C++ and Intel 12.x Fortran builds
#

# Include first to define MKLROOT
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/mkl-12.0.4-options.cmake)

# Set up the basic compilers paths
SET(TRILINOS_TOOLSET_BASE  /opt/gcc-4.5.1/trilinos-toolset)
SET(INTEL_COMPILER_BASE  /opt/intel/Compiler/composerxe-2011.4.191)
# Add rpath for gnu and intel compiler libraries
SET(Trilinos_EXTRA_LINK_FLAGS "-Wl,-rpath,${TRILINOS_TOOLSET_BASE}/lib64" CACHE STRING "")

# Point to path to MKL Fortran modules compatible with GCC 4.5.1 C/C++ calling code
SET(BLAS_INCLUDE_DIRS  ${MKLROOT}/include/intel64/lp64  CACHE PATH  "")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}  CACHE PATH  "")

# To avoid problem with EpetraExt_inout_test failure in optimized code for
# hybrid build
SET(Epetra_ENABLE_Fortran OFF CACHE BOOL "")

# Point to basic CASL-related TPLs related to the GCC C/C++ 4.5.1 compiler
SET(PVMLibraries_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/pvm3/lib/LINUX64 CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/hdf5-1.8.5-patch1/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/hdf5-1.8.5-patch1/include CACHE FILEPATH "")
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "")

# Include last so that above override these cache variables
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-vri-tpls.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-core-enables-disables.cmake)
