INCLUDE(TrilinosCreateClientTemplateHeaders)

# Include Tpetra ETI functions for separate package builds
INCLUDE(${ThyraTpetraAdapters_SOURCE_DIR}/cmake/ETI_functions.cmake)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
APPEND_SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )

# Automatically generate ETI (explicit template instantiation) files
# for Node types that use the Kokkos Devices.
IF(Thyra_ENABLE_EXPLICIT_INSTANTIATION)
  GLOBAL_SET(THYRATPETRAADAPTERS_ETI_CLASSES
    TpetraEuclideanScalarProd
    TpetraLinearOp
    TpetraMultiVector
    TpetraThyraWrappers
    TpetraVector
    TpetraVectorSpace
    )

  # SLGN templates
  TPETRA_PROCESS_ETI_TEMPLATES_SLGN(SRCS
    "ThyraTpetraAdapters_ETI_SC_LO_GO_NT.tmpl"
    "${THYRATPETRAADAPTERS_ETI_CLASSES}"
    "${ThyraTpetraAdapters_ETI_SCALARS}" "${ThyraTpetraAdapters_ETI_LORDS}" "${ThyraTpetraAdapters_ETI_GORDS}" "${ThyraTpetraAdapters_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${SRCS})

ENDIF()


TRIBITS_ADD_LIBRARY(
  thyratpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
