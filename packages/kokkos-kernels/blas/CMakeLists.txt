list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/impl)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/blas/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/tpls)

if(NOT KokkosKernels_ENABLE_COMPONENT_BATCHED)
  message("blas enabled and batched not enabled, we need to include some headers manually!")
  list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched)
  list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
  list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)
endif()

# Adding unit-tests
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/blas)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/blas)

#######################
#                     #
# Logic for BLAS TPLs #
#                     #
#######################

#Include BLAS, Blas host wrapper
if(KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  #Do NOT add this to include path
  append_glob(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tpls/KokkosBlas_Host_tpl.cpp)
endif()

# Include host blas TPL source file
if(KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL OR KOKKOSKERNELS_ENABLE_TPL_ARMPL OR KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)
  list(APPEND SOURCES blas/tpls/KokkosBlas_Host_tpl.cpp)
endif()

# Include cuda blas TPL source file
if(KOKKOSKERNELS_ENABLE_TPL_CUBLAS OR KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  list(APPEND SOURCES blas/tpls/KokkosBlas_Cuda_tpl.cpp)
endif()

# Include rocm blas TPL source file
if(KOKKOSKERNELS_ENABLE_TPL_ROCBLAS)
  list(APPEND SOURCES blas/tpls/KokkosBlas_Rocm_tpl.cpp)
endif()

##################
#                #
# ETI generation #
#                #
##################

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros

macro(gen_blas_eti FUNC_NAME SUBDIR)
  kokkoskernels_generate_eti(${FUNC_NAME} ${SUBDIR}
    COMPONENTS  blas
    HEADER_LIST ETI_HEADERS
    SOURCE_LIST SOURCES
    TYPE_LISTS  FLOATS LAYOUTS DEVICES)
endmacro()

gen_blas_eti(Blas1_abs            abs)
gen_blas_eti(Blas1_abs_mv         abs)
gen_blas_eti(Blas1_scal           scal)
gen_blas_eti(Blas1_scal_mv        scal)
gen_blas_eti(Blas1_dot            dot)
gen_blas_eti(Blas1_dot_mv         dot)
gen_blas_eti(Blas1_axpby          axpby)
gen_blas_eti(Blas1_axpby_mv       axpby)
gen_blas_eti(Blas1_update         update)
gen_blas_eti(Blas1_update_mv      update)
gen_blas_eti(Blas1_sum            sum)
gen_blas_eti(Blas1_sum_mv         sum)
gen_blas_eti(Blas1_nrm1           nrm1)
gen_blas_eti(Blas1_nrm1_mv        nrm1)
gen_blas_eti(Blas1_nrm2w          nrm2w)
gen_blas_eti(Blas1_nrm2w_mv       nrm2w)
gen_blas_eti(Blas1_nrminf         nrminf)
gen_blas_eti(Blas1_nrminf_mv      nrminf)
gen_blas_eti(Blas1_iamax          iamax)
gen_blas_eti(Blas1_iamax_mv       iamax)
gen_blas_eti(Blas1_nrm2           nrm2)
gen_blas_eti(Blas1_nrm2_mv        nrm2)
gen_blas_eti(Blas1_mult           mult)
gen_blas_eti(Blas1_mult_mv        mult)
gen_blas_eti(Blas1_reciprocal     reciprocal)
gen_blas_eti(Blas1_reciprocal_mv  reciprocal)
gen_blas_eti(Blas1_rot            rot)
gen_blas_eti(Blas1_rotg           rotg)
gen_blas_eti(Blas1_swap           swap)
gen_blas_eti(Blas2_gemv           gemv)
gen_blas_eti(Blas2_ger            ger)
gen_blas_eti(Blas2_syr            syr)
gen_blas_eti(Blas2_syr2           syr2)
gen_blas_eti(Blas3_gemm           gemm)
gen_blas_eti(Blas3_trsm           trsm)
gen_blas_eti(Blas3_trmm           trmm)

# These last two do REAL_FLOATS

kokkoskernels_generate_eti(Blas1_rotm rotm
  COMPONENTS  blas
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  REAL_FLOATS LAYOUTS DEVICES
)

kokkoskernels_generate_eti(Blas1_rotmg rotmg
  COMPONENTS  blas
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  REAL_FLOATS LAYOUTS DEVICES
)
