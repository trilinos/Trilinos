
SET(INTEL_11_1_064_ROOT /opt/intel/Compiler/11.1/064)
SET(MKLROOT ${INTEL_11_1_064_ROOT}/mkl)
SET(MKL_IFORT_MODULE_PATH ${MKLROOT}/include/em64t/lp64)
SET(MKLLIB "${MKLROOT}/lib/em64t")

SET(BLAS_LIBRARY_NAMES  "mkl_intel_lp64;mkl_blas95_lp64;mkl_core;mkl_sequential"
  CACHE STRING  "Sequential 64-bit BLAS with 32-bit integers")
SET(BLAS_LIBRARY_DIRS  "${MKLLIB}"  CACHE FILEPATH "")
SET(LAPACK_LIBRARY_NAMES "mkl_lapack95_lp64"
  CACHE STRING  "64-bit LAPACK with 32-bit integers")
SET(LAPACK_LIBRARY_DIRS  "${MKLLIB}"  CACHE FILEPATH "")
