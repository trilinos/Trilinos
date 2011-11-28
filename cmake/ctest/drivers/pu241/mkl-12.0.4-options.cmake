
SET(INTEL_12_0_4_ROOT /opt/intel/Compiler/composerxe-2011.4.191)
SET(MKLROOT ${INTEL_12_0_4_ROOT}/mkl)
SET(MKL_IFORT_MODULE_PATH ${MKLROOT}/include/intel64/lp64)
SET(MKL_GCC451_MODULE_PATH ${MKLROOT}/include/gcc-4.5.1)
SET(MKLLIB "${MKLROOT}/lib/intel64")

SET(BLAS_LIBRARY_NAMES  "mkl_intel_lp64;mkl_blas95_lp64;mkl_core;mkl_sequential"
  CACHE STRING  "Sequential 64-bit BLAS with 32-bit integers")
SET(BLAS_LIBRARY_DIRS  "${MKLLIB}"  CACHE FILEPATH "")
SET(LAPACK_LIBRARY_NAMES "mkl_lapack95_lp64"
  CACHE STRING  "64-bit LAPACK with 32-bit integers")
SET(LAPACK_LIBRARY_DIRS  "${MKLLIB}"  CACHE FILEPATH "")
