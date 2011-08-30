INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/mkl-12.0.4-options.cmake)
SET(BLAS_INCLUDE_DIRS   ${MKLROOT}/include/intel64/lp64                     CACHE PATH  "Path to MKL BLAS Fortran modules compatible with Intel fortran")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}                                CACHE PATH  "Path to MKL LAPACK Fortran modules compatible with Intel fortran")

INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/boost-1.46.1-options.cmake)
INCLUDE(${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/casl-vri-tpls.cmake)

SET(PVM_LIBRARY_DIRS /opt/intel-11.1.064/tpls/pvm3/lib/LINUX64                    CACHE FILEPATH "")
SET(PVM_INCLUDE_DIRS /opt/intel-11.1.064/tpls/pvm3/include                        CACHE FILEPATH "")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;hdf5_cpp"                                    CACHE STRING   "")
SET(HDF5_LIBRARY_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/lib              CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS /opt/intel-11.1.064/tpls/hdf5-1.8.5-patch1/include          CACHE FILEPATH "")
SET(EpetraExt_ENABLE_HDF5 OFF CACHE BOOL "")

SET(INTEL_BIN /opt/intel/Compiler/composerxe-2011.4.191/bin/intel64)

SET(CMAKE_SKIP_RPATH ON BOOL "")
#SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
SET(TPL_ENABLE_BinUtils ON CACHE BOOL "")
SET(CMAKE_C_COMPILER ${INTEL_BIN}/icc CACHE FILEPATH "")
SET(CMAKE_CXX_COMPILER ${INTEL_BIN}/icpc CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER ${INTEL_BIN}/ifort CACHE FILEPATH "")
