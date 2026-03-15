function(kokkoskernels_feature_depends_on_tpls FEATURE)
  if(KOKKOSKERNELS_ENABLE_${FEATURE})
    foreach(TPL ${ARGN})
      if(NOT KOKKOSKERNELS_ENABLE_TPL_${TPL})
        message(SEND_ERROR
          "Feature ${FEATURE} requires TPL support for ${TPL}. Must build with -DKokkosKernels_ENABLE_TPL_${TPL}:BOOL=ON and potentially -D${TPL}_ROOT=<INSTALL> to the desired package location")
      endif()
    endforeach()
  endif()
endfunction()

kokkoskernels_add_option("ENABLE_SUPERNODAL_SPTRSV" ON BOOL
  "Whether to build supernodal SPTRSV support")

if(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV AND NOT KOKKOSKERNELS_INST_LAYOUTLEFT)
  message(WARNING "Disabling SUPERNODAL_SPTRSV - this capability is only supported with LayoutLeft")
  set(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV OFF CACHE BOOL
     "Disabling SUPERNODAL_SPTRSV - this capability is only supported with LayoutLeft"
     FORCE)
endif()

# ==================================================================
# Fortran Complex BLAS
# ==================================================================

if(KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  include(CheckHostBlasReturnComplex.cmake)
  check_host_blas_return_complex(KOKKOSKERNELS_TPL_BLAS_RETURN_COMPLEX)
endif()

# ==================================================================
# Lapack requirements
# ==================================================================

if(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_ROCBLAS AND NOT KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
  message(FATAL_ERROR
    "rocSOLVER requires rocBLAS and rocSPARSE, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_ROCBLAS:BOOL=ON and KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE:BOOL=ON.")
elseif(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
  message(FATAL_ERROR
    "rocSOLVER requires rocSPARSE, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE:BOOL=ON.")
elseif(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
  message(FATAL_ERROR
    "rocSOLVER requires rocBLAS, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_ROCBLAS:BOOL=ON.")
endif()

# TPL_ENABLE_CUDA default enables CUBLAS and CUSOLVER in Trilinos, but not CUSPARSE. CUSPARSE is a required TPL for CUSOLVER support in KokkosKernels.
if(KOKKOSKERNELS_HAS_TRILINOS AND TPL_ENABLE_CUDA)
  # Checks disable CUSOLVER in KokkosKernels if TPL dependency requirements are not met. This is a compatibility workaround to allow existing configuration options for Trilinos to continue working.
  if(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUBLAS AND NOT KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    message(WARNING
      "cuSOLVER requires cuBLAS and cuSPARSE, disabling cuSOLVER. To use cuSOLVER, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUBLAS:BOOL=ON and KOKKOSKERNELS_ENABLE_TPL_CUSPARSE:BOOL=ON to use.")
    set(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER OFF CACHE BOOL
       "Disabling KOKKOSKERNELS_ENABLE_TPL_CUSOLVER - this capability requires both CUBLAS and CUSPARSE TPLs"
       FORCE)
  elseif(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    message(WARNING
      "cuSOLVER requires cuSPARSE, disabling cuSOLVER. To use cuSOLVER, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUSPARSE:BOOL=ON to use.")
    set(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER OFF CACHE BOOL
        "Disabling KOKKOSKERNELS_ENABLE_TPL_CUSOLVER - this capability requires both CUBLAS and CUSPARSE TPLs"
        FORCE)
  elseif(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    message(WARNING
      "cuSOLVER requires cuBLAS, disabling cuSOLVER. To use cuSOLVER, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUBLAS:BOOL=ON to use.")
    set(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER OFF CACHE BOOL
       "Disabling KOKKOSKERNELS_ENABLE_TPL_CUSOLVER - this capability requires both CUBLAS and CUSPARSE TPLs"
       FORCE)
  endif()
else()
  if(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUBLAS AND NOT KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    message(FATAL_ERROR
      "cuSOLVER requires cuBLAS and cuSPARSE, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUBLAS:BOOL=ON and KOKKOSKERNELS_ENABLE_TPL_CUSPARSE:BOOL=ON.")
  elseif(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
    message(FATAL_ERROR
      "cuSOLVER requires cuSPARSE, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUSPARSE:BOOL=ON.")
  elseif(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER AND NOT KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    message(FATAL_ERROR
      "cuSOLVER requires cuBLAS, please reconfigure with KOKKOSKERNELS_ENABLE_TPL_CUBLAS:BOOL=ON.")
  endif()
endif()
