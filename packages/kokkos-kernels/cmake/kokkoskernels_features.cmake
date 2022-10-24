function(kokkoskernels_feature_depends_on_tpls FEATURE)
  IF (KOKKOSKERNELS_ENABLE_${FEATURE})
    FOREACH(TPL ${ARGN})
      IF (NOT KOKKOSKERNELS_ENABLE_TPL_${TPL})
        MESSAGE(SEND_ERROR "Feature ${FEATURE} requires TPL support for ${TPL}. Must build with -DKokkosKernels_ENABLE_TPL_${TPL}:BOOL=ON and potentially -D${TPL}_ROOT=<INSTALL> to the desired package location")
      ENDIF()
    ENDFOREACH()
  ENDIF()
endfunction()

KOKKOSKERNELS_ADD_OPTION(
  ENABLE_SUPERNODAL_SPTRSV
  ON
  BOOL
  "Whether to build supernodal SPTRSV support")

IF (KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV AND NOT KOKKOSKERNELS_INST_LAYOUTLEFT)
  MESSAGE(WARNING "Disabling SUPERNODAL_SPTRSV - this capability is only supported with LayoutLeft")
  SET(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV OFF CACHE BOOL "Disabling SUPERNODAL_SPTRSV - this capability is only supported with LayoutLeft" FORCE)
ENDIF()

# ==================================================================
# Fortran Complex BLAS
# ==================================================================

IF (KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  INCLUDE(CheckHostBlasReturnComplex.cmake)
  CHECK_HOST_BLAS_RETURN_COMPLEX(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
ENDIF()
