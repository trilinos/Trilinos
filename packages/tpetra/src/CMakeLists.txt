INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

# Add Kokkos refactor code
IF (Tpetra_Have_Kokkos_Refactor)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_refactor)
  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})
ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (KokkosClassic_ENABLE_Thrust)
  FILE(GLOB CUDA_SOURCES *.cu)
  PACKAGE_ADD_CUDA_LIBRARY(tpetra_cuda SOURCES ${CUDA_SOURCES})
ENDIF()

TRIBITS_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )