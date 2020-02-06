
KOKKOSKERNELS_CONFIGURE_FILE(KokkosKernels_config.h)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

#I don't want the complexity of sometimes interface/sometimes not
#Just minimally force an empty dummy file
SET(SOURCES dummy.cpp)

#-----------------------------------------------------------------------------

LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/impl)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/impl/tpls)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/blas/impl)

#Include Sparse
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sparse/impl)

#Include graph
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/impl)

#Include batched
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched)

#Include Common
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/common)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/impl)


FOREACH(DIR ${KK_INCLUDE_DIRS})
  KOKKOSKERNELS_INCLUDE_DIRECTORIES(${DIR})
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
ENDFOREACH()


#Include BLAS, Blas host wrapper
IF (KOKKOSKERNELS_ENABLE_TPL_BLAS OR KOKKOSKERNELS_ENABLE_TPL_MKL)
  #Do NOT add this to include path
  APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/impl/tpls/KokkosBlas_Host_tpl.cpp)
ENDIF()

SET(IMPL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)
IF (KOKKOSKERNELS_INST_EXECSPACE_SERIAL)
  IF (KOKKOSKERNELS_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_dbl_*_Serial_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_flt_*_Serial_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_dbl_*_Serial_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_flt_*_Serial_*.cpp)
  ENDIF()
ENDIF()

IF (KOKKOSKERNELS_INST_EXECSPACE_THREADS)
  IF (KOKKOSKERNELS_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_dbl_*_Threads_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_flt_*_Threads_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_dbl_*_Threads_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_flt_*_Threads_*.cpp)
  ENDIF()
ENDIF()

IF (KOKKOSKERNELS_INST_EXECSPACE_OPENMP)
  IF (KOKKOSKERNELS_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_dbl_*_OpenMP_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_flt_*_OpenMP_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_dbl_*_OpenMP_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_flt_*_OpenMP_*.cpp)
  ENDIF()
ENDIF()

IF (KOKKOSKERNELS_INST_EXECSPACE_CUDA)
  APPEND_GLOB(SOURCES ${IMPL_DIR}/tpls/*_Cuda_*.cpp)
  IF (KOKKOSKERNELS_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_dbl_*_Cuda_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*inst_flt_*_Cuda_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_dbl_*_Cuda_*.cpp)
  ENDIF()
  IF (KOKKOSKERNELS_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${IMPL_DIR}/generated_specializations_cpp/*/*_cmplx_flt_*_Cuda_*.cpp)
  ENDIF()
ENDIF()

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)
#LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

#-----------------------------------------------------------------------------

KOKKOSKERNELS_ADD_LIBRARY(
  kokkoskernels
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)


IF (KOKKOSKERNELS_HAS_TRILINOS)
  #no linking commands required - tribits does this
ELSE()
  TARGET_LINK_LIBRARIES(kokkoskernels PUBLIC ${KOKKOS_DEPENDENCY_NAME})
  FOREACH(DIR ${KK_INCLUDE_DIRS})
    TARGET_INCLUDE_DIRECTORIES(kokkoskernels PUBLIC $<BUILD_INTERFACE:${DIR}>)
  ENDFOREACH()
  TARGET_INCLUDE_DIRECTORIES(kokkoskernels PUBLIC $<INSTALL_INTERFACE:include>)
ENDIF()

KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC BLAS)
KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC LAPACK)
KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC MKL)
KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CUBLAS)
KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC CUSPARSE)
# Not yet here KOKKOSKERNELS_LINK_TPL(kokkoskernels PUBLIC MAGMA)

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/impl/generated_specializations_hpp)
APPEND_GLOB(HEADERS_GENERATED ${CMAKE_CURRENT_SOURCE_DIR}/impl/generated_specializations_hpp/*.hpp)
INSTALL(FILES ${HEADERS_GENERATED} DESTINATION ${KOKKOSKERNELS_HEADER_DIR}/generated_specializations_hpp/)

#
# Make a trivial change here if you add a file either to this
# directory, or to the impl/ subdirectory.  Keep making changes
# whenever you add or remove a .cpp or .hpp file, and everything will
# be OK.  Here is another such change, and here is yet another.
#

