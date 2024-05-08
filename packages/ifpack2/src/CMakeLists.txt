INCLUDE(TrilinosCreateClientTemplateHeaders)

# Parse the Hypre headers
IF (${PACKAGE_NAME}_ENABLE_HYPRE)
  IF (PYTHON_EXECUTABLE)
    EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../utils/parseHypre.py ${HYPRE_INCLUDE_DIRS} Ifpack2_HypreParameterMap.hpp
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


FUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE_LGN ETI_CLASSES TEMPLATE_FILE SOURCES_LIST IFPACK2_ETI_LORDS IFPACK2_ETI_GORDS IFPACK2_ETI_NODES)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(CLASS_MANGLED ${CLASS})
    string(TOUPPER "${CLASS_MANGLED}" UPPER_CASE_CLASS)
    TPETRA_PROCESS_ALL_LGN_TEMPLATES(TMP_OUTPUT_FILES ${TEMPLATE_FILE} ${CLASS_MANGLED} ${UPPER_CASE_CLASS} "${IFPACK2_ETI_LORDS}" "${IFPACK2_ETI_GORDS}" "${IFPACK2_ETI_NODES}")
    LIST(APPEND SRCS ${TMP_OUTPUT_FILES})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE_LGN)

FUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE_SLGN ETI_CLASSES TEMPLATE_FILE SOURCES_LIST IFPACK2_ETI_SCALARS IFPACK2_ETI_LORDS IFPACK2_ETI_GORDS IFPACK2_ETI_NODES)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(CLASS_MANGLED ${CLASS})
    string(TOUPPER "${CLASS_MANGLED}" UPPER_CASE_CLASS)
    TPETRA_PROCESS_ALL_SLGN_TEMPLATES(TMP_OUTPUT_FILES ${TEMPLATE_FILE} ${CLASS_MANGLED} ${UPPER_CASE_CLASS} "${IFPACK2_ETI_SCALARS}" "${IFPACK2_ETI_LORDS}" "${IFPACK2_ETI_GORDS}" "${IFPACK2_ETI_NODES}" FALSE)
    LIST(APPEND SRCS ${TMP_OUTPUT_FILES})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE_SLGN)

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
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# SupportGraph requires some TPLs and options to be enabled.
IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Lemon AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Cholmod)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/supportgraph)
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
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h )

# Automatically generate ETI (explicit template instanatiation) files
# for Node types that use the (new) Kokkos Devices.
SET(IFPACK2_CPP_SOURCES "")
IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
  # Set the list of Ifpack2 classes templated on <Scalar, LO, GO,
  # Node> for which we want to do ETI using this system.  These
  # classes usually operate on sparse matrices (instances of
  # Tpetra::CrsMatrix or Tpetra::RowMatrix, which also take these
  # template parameters).
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

  GLOBAL_SET(IFPACK2_ETI_LINEAR_SOLVER_FACTORY_CLASSES
    Details::LinearSolverFactory
  )

  SET(TEMPLATE_FILE "Ifpack2_Details_LinearSolverFactory_ETI_SC_LO_GO_NT.tmpl")
  IFPACK2_PROCESS_ETI_TEMPLATE_SLGN(${IFPACK2_ETI_LINEAR_SOLVER_FACTORY_CLASSES} ${TEMPLATE_FILE} IFPACK2_FACTORY_OUTPUT_FILES "${Ifpack2_ETI_SCALARS}" "${Ifpack2_ETI_LORDS}" "${Ifpack2_ETI_GORDS}" "${Ifpack2_ETI_NODES}")
  LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_FACTORY_OUTPUT_FILES})

  SET(TEMPLATE_FILE "Ifpack2_ETI_LO_GO_NT.tmpl")
  IFPACK2_PROCESS_ETI_TEMPLATE_LGN("${IFPACK2_ETI_LO_GO_CLASSES}" ${TEMPLATE_FILE} IFPACK2_SRCS "${Ifpack2_ETI_LORDS}" "${Ifpack2_ETI_GORDS}" "${Ifpack2_ETI_NODES}")
  LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})

  SET(TEMPLATE_FILE "Ifpack2_ETI_SC_LO_GO_NT.tmpl")
  IFPACK2_PROCESS_ETI_TEMPLATE_SLGN("${IFPACK2_ETI_CLASSES}" ${TEMPLATE_FILE} IFPACK2_SRCS "${Ifpack2_ETI_SCALARS}" "${Ifpack2_ETI_LORDS}" "${Ifpack2_ETI_GORDS}" "${Ifpack2_ETI_NODES}")
  LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})

ENDIF()

#MESSAGE(DEBUG " *** IFPACK2_CPP_SOURCES = ${IFPACK2_CPP_SOURCES}")

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${IFPACK2_CPP_SOURCES}
  )

#
# Make a trivial change here if you want CMake to run, due to changes
# you make to files in Ifpack2.  Here is another such change.
# Behold, I make another change, and another, and another.
#
