
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(HEADERS ${DIR}/generated_specializations_hpp/*.hpp)

IF (KokkosKernels_INST_EXECSPACE_SERIAL)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_double_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_float_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Serial_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Serial_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_PTHREAD)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_double_*_Pthread_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_float_*_Pthread_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Pthread_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Pthread_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_OPENMP)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_double_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_float_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_OpenMP_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_OpenMP_*.cpp)
  ENDIF()
ENDIF()

IF (KokkosKernels_INST_EXECSPACE_CUDA)
  IF (KokkosKernels_INST_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_double_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*specialization_float_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_DOUBLE)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_double_*_Cuda_*.cpp)
  ENDIF()
  IF (KokkosKernels_INST_COMPLEX_FLOAT)
    APPEND_GLOB(SOURCES ${DIR}/generated_specializations_cpp/*/*_complex_float_*_Cuda_*.cpp)
  ENDIF()
ENDIF()

LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

#-----------------------------------------------------------------------------

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stage/graph)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stage/graph/impl)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stage/graph/utils)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

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

