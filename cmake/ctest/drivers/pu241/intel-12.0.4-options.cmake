#
# Core options for Intel 12.0.4 builds
#

# TPLs that are compiled with the Intel C/C++ compiler
SET(PVMLibraries_LIBRARY_DIRS /opt/intel-11.1.064/tpls/pvm3/lib/LINUX64 CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS /opt/intel-11.1.064/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/include CACHE FILEPATH "")

# Define the core compilers
SET(INTEL_BIN /opt/intel/Compiler/composerxe-2011.4.191/bin/intel64)
SET(CMAKE_C_COMPILER ${INTEL_BIN}/icc CACHE FILEPATH "")
SET(CMAKE_CXX_COMPILER ${INTEL_BIN}/icpc CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER ${INTEL_BIN}/ifort CACHE FILEPATH "")

SET(CMAKE_SKIP_RPATH ON BOOL "")
#SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")

# Down at the bottom to allow the above cache varibles to override
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/mkl-12.0.4-options.cmake)

# Must be set below above include because it defines MPKROOT
SET(BLAS_INCLUDE_DIRS ${MKLROOT}/include/intel64/lp64 CACHE PATH "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

# Including these last to allow override above
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-vri-tpls.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-core-enables-disables.cmake)
