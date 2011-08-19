SET(INTEL_12_0_4_ROOT /opt/intel/Compiler/composerxe-2011.4.191)

SET(MKLROOT ${INTEL_12_0_4_ROOT}/mkl)

SET(BLAS_LIBRARY_NAMES   "mkl_intel_lp64;mkl_blas95_lp64;mkl_core;mkl_sequential"   CACHE STRING   "")
SET(BLAS_LIBRARY_DIRS    "${MKLROOT}/lib/intel64"                                   CACHE FILEPATH "")
SET(LAPACK_LIBRARY_NAMES "mkl_lapack95_lp64"                                        CACHE STRING   "")
SET(LAPACK_LIBRARY_DIRS  "${MKLROOT}/lib/intel64"                                   CACHE FILEPATH "")

# hijack these to get modules paths into the include list; needed by CASLRAVE/ANC
SET(BLAS_INCLUDE_DIRS   ${MKLROOT}/include/intel64/lp64                     CACHE PATH     "")
SET(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS}                                CACHE PATH     "")
