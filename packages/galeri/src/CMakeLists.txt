INCLUDE(TrilinosCreateClientTemplateHeaders)

# Include Tpetra ETI functions for separate package builds
INCLUDE(${Galeri_SOURCE_DIR}/cmake/ETI_functions.cmake)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../headers)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

APPEND_SET(HEADERS
  ../headers/Galeri_ConfigDefs.h
  ../headers/Galeri_Exception.h

  Galeri_XpetraProblemFactory.hpp
  Galeri_XpetraProblemFactory_Helmholtz.hpp
  Galeri_XpetraMatrixTypes.hpp
  Galeri_XpetraMatrixTypes_Helmholtz.hpp
  Galeri_StencilProblems_Helmholtz.hpp
  Galeri_Problem_Helmholtz.hpp
  Galeri_HelmholtzFEM2DProblem.hpp
  Galeri_HelmholtzFEM3DProblem.hpp
  Galeri_VelocityModel.hpp

  Galeri_VectorTraits.hpp
  Galeri_MatrixTraits.hpp
  Galeri_MapTraits.hpp
  Galeri_MultiVectorTraits.hpp

  Galeri_XpetraUtils.hpp
  Galeri_Problem.hpp
  )

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*_decl.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*_def.hpp)

APPEND_SET(SOURCES
  mac.cpp
  )

# Automatically generate ETI (explicit template instantiation) files
# for Node types that use the Kokkos Devices.
IF(Galeri_ENABLE_EXPLICIT_INSTANTIATION)
  # Set the list of Galeri classes templated on <Scalar, LO, GO,
  # Node> for which we want to do ETI using this system.  These
  # classes usually operate on sparse matrices (instances of
  # Tpetra::CrsMatrix or Tpetra::RowMatrix, which also take these
  # template parameters).

  # GALERI_ETI_CLASSES and GALERI_ETI_LO_GO_CLASSES are global
  # variables so that Stokhos can do ETI.

  GLOBAL_SET(GALERI_ETI_CLASSES
    Elasticity2DProblem
    Elasticity3DProblem
    StencilProblems
  )

  GLOBAL_SET(GALERI_ETI_LO_GO_CLASSES
    XpetraCartesian
    XpetraMaps
    )

  GLOBAL_SET(GALERI_ETI_GO_CLASSES
    XpetraParameters
    )

  # G templates
  TPETRA_PROCESS_ETI_TEMPLATES_G(SRCS
    "Galeri_ETI_GO.tmpl"
    "${GALERI_ETI_GO_CLASSES}"
    "${Galeri_ETI_GORDS}")
  LIST(APPEND SOURCES ${SRCS})

  # LGN templates
  TPETRA_PROCESS_ETI_TEMPLATES_LGN(SRCS
    "Galeri_ETI_LO_GO_NT.tmpl"
    "${GALERI_ETI_LO_GO_CLASSES}"
    "${Galeri_ETI_LORDS}" "${Galeri_ETI_GORDS}" "${Galeri_ETI_NODES}")
  LIST(APPEND SOURCES ${SRCS})

  # SLGN templates
  TPETRA_PROCESS_ETI_TEMPLATES_SLGN(SRCS
    "Galeri_ETI_SC_LO_GO_NT.tmpl"
    "${GALERI_ETI_CLASSES}"
    "${Galeri_ETI_SCALARS}" "${Galeri_ETI_LORDS}" "${Galeri_ETI_GORDS}" "${Galeri_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${SRCS})

ENDIF()

TRIBITS_ADD_LIBRARY(
  galeri-xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
