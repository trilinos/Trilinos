INCLUDE(TrilinosCreateClientTemplateHeaders)

# Parse the Hypre headers
IF (${PACKAGE_NAME}_ENABLE_HYPRE)
  IF (Python3_EXECUTABLE)
    EXECUTE_PROCESS(
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../utils/parseHypre.py ${HYPRE_INCLUDE_DIRS} Ifpack2_HypreParameterMap.hpp
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      RESULT_VARIABLE RETURN_VALUE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  ELSE ()
    MESSAGE(FATAL_ERROR "Python needed to parse Hypre headers")
  ENDIF ()
  IF (RETURN_VALUE EQUAL 0)
    MESSAGE("-- Parsed Hypre headers")
  ELSE ()
    MESSAGE(FATAL_ERROR "Failed to parse Hypre headers")
  ENDIF ()
ENDIF ()

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

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})

# SupportGraph requires some TPLs and options to be enabled.
IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Lemon AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Cholmod)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/supportgraph)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/supportgraph)
  SET(HEADERS ${HEADERS}
    supportgraph/Ifpack2_SupportGraph_decl.hpp
    supportgraph/Ifpack2_SupportGraph_def.hpp
    )
  SET(SOURCES ${SOURCES}
    supportgraph/Ifpack2_SupportGraph.cpp
    )
ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
APPEND_SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )
# APPEND_SET(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_ETIHelperMacros.h )

# Automatically generate ETI (explicit template instantiation) files
# for Node types that use the Kokkos Devices.
IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
  # Set the list of Ifpack2 classes templated on <Scalar, LO, GO,
  # Node> for which we want to do ETI using this system.  These
  # classes usually operate on sparse matrices (instances of
  # Tpetra::CrsMatrix or Tpetra::RowMatrix, which also take these
  # template parameters).

  # IFPACK2_ETI_CLASSES and IFPACK2_ETI_LO_GO_CLASSES are global
  # variables so that Stokhos can do ETI.

  GLOBAL_SET(IFPACK2_ETI_CLASSES
    AdditiveSchwarz
    BlockRelaxation
    BorderedOperator
    Chebyshev
    Container
    BandedContainer
    DatabaseSchwarz
    DenseContainer
    Diagonal
    DiagonalFilter
    DropFilter
    Factory
    Hiptmair
    IdentitySolver
    ILUT
    LinePartitioner
    LocalFilter
    LocalSparseTriangularSolver
    MDF
    OverlappingRowMatrix
    Relaxation
    ReorderFilter
    RILUK
    SingletonFilter
    SparseContainer
    SparsityFilter
    ContainerFactory
    TriDiContainer
    Details::AdditiveSchwarzFilter
    Details::Chebyshev
    Details::ChebyshevKernel
    Details::DenseSolver
    Details::Factory
    Details::InverseDiagonalKernel
    Details::LinearSolver
    Details::OneLevelFactory
    Details::ScaledDampedResidual
    Details::TriDiSolver
    Experimental::RBILUK
    Details::LinearSolverFactory
    )

  IF(Ifpack2_ENABLE_Amesos2)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Details::Amesos2Wrapper
      )
  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_Experimental_KokkosKernels_Features)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES BlockTriDiContainer)
  ENDIF()

  IF(Ifpack2_ENABLE_ShyLU_NodeFastILU)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Details::FastILU_Base
      Details::Filu
      Details::Fic
      Details::Fildl
    )
  ENDIF()

  IF(Ifpack2_ENABLE_HYPRE)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Hypre
      )
  ENDIF()

  # Set the list of Ifpack2 classes templated on <LO, GO, Node> for
  # which we want to do ETI using this system.  These classes usually
  # operate on sparse graphs (instances of Tpetra::CrsGraph or
  # Tpetra::RowGraph, which also take these template parameters).
  GLOBAL_SET(IFPACK2_ETI_LO_GO_CLASSES
    LinearPartitioner
    OverlappingPartitioner
    Details::OverlappingRowGraph
    )

  # LGN templates
  TPETRA_PROCESS_ETI_TEMPLATES_LGN(SRCS
    "Ifpack2_ETI_LO_GO_NT.tmpl"
    "${IFPACK2_ETI_LO_GO_CLASSES}"
    "${Ifpack2_ETI_LORDS}" "${Ifpack2_ETI_GORDS}" "${Ifpack2_ETI_NODES}")
  LIST(APPEND SOURCES ${SRCS})

  # SLGN templates
  TPETRA_PROCESS_ETI_TEMPLATES_SLGN(SRCS
    "Ifpack2_ETI_SC_LO_GO_NT.tmpl"
    "${IFPACK2_ETI_CLASSES}"
    "${Ifpack2_ETI_SCALARS}" "${Ifpack2_ETI_LORDS}" "${Ifpack2_ETI_GORDS}" "${Ifpack2_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${SRCS})

ENDIF()

#MESSAGE(DEBUG " *** IFPACK2_CPP_SOURCES = ${IFPACK2_CPP_SOURCES}")

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# Make a trivial change here if you want CMake to run, due to changes
# you make to files in Ifpack2.  Here is another such change.
# Behold, I make another change, and another, and another.
#
