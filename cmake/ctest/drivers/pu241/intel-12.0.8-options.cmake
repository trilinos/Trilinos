#
# Core options for Intel 12.0.8 builds
#

# Define the core compilers
SET(IFORT_VERSION "12.0.8")
SET(INTEL_LIB /opt/intel-12.0.8/lib/intel64)
SET(INTEL_BIN /opt/intel-12.0.8/bin)
# Add rpath for compiler libraries
# SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS "-Wl,-rpath,${INTEL_LIB}" CACHE STRING "")
SET(CMAKE_SKIP_RPATH ON BOOL "")

# this compiler supports BinUtils
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")

# Include MKL and TBB; these should match version of Intel compilers being used
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/tbb-12.0.8-options.cmake)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/mkl-12.0.8-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${MKL_IFORT_MODULE_PATH} CACHE PATH "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

# Point to basic CASL-related TPLs compatible with the Intel C/C++ compiler
SET(PVMLibraries_LIBRARY_DIRS /opt/intel-12.0.8/tpls/pvm3/lib     CACHE FILEPATH "")
SET(PVMHeaders_INCLUDE_DIRS   /opt/intel-12.0.8/tpls/pvm3/include CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp" CACHE STRING "")
SET(HDF5_LIBRARY_DIRS /opt/intel-12.0.8/tpls/hdf5-1.8.5-patch1/lib     CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/intel-12.0.8/tpls/hdf5-1.8.5-patch1/include CACHE FILEPATH "")

# Including these last to allow override above
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/casl-core-enables-disables.cmake)
