
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(HEADERS ${DIR}/generated_specializations_hpp/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl/tpls)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#Include BLAS

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/blas)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/blas/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#Include Sparse

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sparse)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sparse/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#Include graph

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/graph)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/graph/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#Include batched

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/batched)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/batched/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#Include Common

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/common)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)

IF (KokkosKernels_INST_EXECSPACE_SERIAL)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_double_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_float_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Serial_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_THREADS)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_double_*_Threads_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_float_*_Threads_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Threads_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Threads_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_OPENMP)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_double_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_float_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_OpenMP_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_CUDA)
  APPEND_GLOB(SOURCES ${DIR}/tpls/*_Cuda_*.cpp)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_double_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*inst_float_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Cuda_*.cpp)
  ENDIF()
ENDIF()

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)
#LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
  kokkoskernels
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl/generated_specializations_hpp)
APPEND_GLOB(HEADERS_GENERATED  ${DIR}/*.hpp)
INSTALL(FILES ${HEADERS_GENERATED} DESTINATION ${TRILINOS_INCDIR}/generated_specializations_hpp/)
#
# Make a trivial change here if you add a file either to this
# directory, or to the impl/ subdirectory.  Keep making changes
# whenever you add or remove a .cpp or .hpp file, and everything will
# be OK.  Here is another such change, and here is yet another.
#

