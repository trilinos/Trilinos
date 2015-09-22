#
# Core options for Intel 11.1.064 builds
#

# Define the core compilers
SET(IFORT_VERSION "11.1.064")
SET(INTEL_LIB /opt/intel/Compiler/11.1/064/lib/intel64)
SET(INTEL_BIN /opt/intel/Compiler/11.1/064/bin/intel64)
# Add rpath for compiler libraries
# SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS "-Wl,-rpath,${INTEL_LIB}" CACHE STRING "")
SET(CMAKE_SKIP_RPATH ON BOOL "")
# Intel 11 not compatible with GCC 4.5.1; in case it is loaded, lock it out
SET(CMAKE_C_FLAGS   "-gcc-name=/usr/bin/gcc" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-gxx-name=/usr/bin/g++" CACHE STRING "")

# this compiler supports BinUtils
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")

# Include MKL and TBB; these should match version of Intel compilers being used
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-11.1.064-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-11.1.064-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

# Point to basic CASL-related TPLs compatible with the Intel C/C++ compiler
SET(PVMLibraries_LIBRARY_DIRS /opt/intel-11.1.064/tpls/pvm3/lib/LINUX64 CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS /opt/intel-11.1.064/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/lib CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/include CACHE FILEPATH "")

# Including these last to allow override above
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
