INCLUDE(TrilinosCreateClientTemplateHeaders)

# Function to generate ETI files from a template and list of class names
FUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(IFPACK2_PROCESS_ETI_TEMPLATE)

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

IF (${PACKAGE_NAME}_ENABLE_Experimental)
IF (${PACKAGE_NAME}_ENABLE_Lemon)
IF (${PACKAGE_NAME}_ENABLE_Amesos2)
IF (${PACKAGE_NAME}_ENABLE_Cholmod)

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/supportgraph)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/supportgraph)

SET(HEADERS ${HEADERS}
supportgraph/Ifpack2_SupportGraph_decl.hpp
supportgraph/Ifpack2_SupportGraph_def.hpp
)
SET(SOURCES ${SOURCES}
supportgraph/Ifpack2_SupportGraph.cpp
)

ENDIF()
ENDIF()
ENDIF()
ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )

# Auto generate ETI files for new Tpetra nodes
SET(IFPACK2_CPP_SOURCES "")
SET(IFPACK2_CU_SOURCES "")
IF(Ifpack2_ENABLE_EXPLICIT_INSTANTIATION)
  GLOBAL_SET(IFPACK2_ETI_CLASSES
    AdditiveSchwarz
    BlockRelaxation
    BorderedOperator
    Chebyshev
    DenseContainer
    Diagonal
    DiagonalFilter
    Factory
    Hiptmair
    IdentitySolver
    ILUT
    Krylov
    LocalFilter
    OverlappingRowMatrix
    Relaxation
    ReorderFilter
    RILUK
    SingletonFilter
    SparseContainer
    SparsityFilter
    Details::Chebyshev
    Details::DenseSolver
    Details::OneLevelFactory
    )
  IF(Ifpack2_ENABLE_Amesos2)
    APPEND_GLOBAL_SET(IFPACK2_ETI_CLASSES
      Details::Amesos2Wrapper
      )
  ENDIF()
  GLOBAL_SET(IFPACK2_ETI_LO_GO_CLASSES
    LinearPartitioner
    OverlappingPartitioner
    Details::OverlappingRowGraph
    )
  IF (Kokkos_ENABLE_Pthread)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Threads.tmpl
      "Ifpack2_CLASS_FILE_NAME_Threads.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Threads.tmpl
      "Ifpack2_CLASS_FILE_NAME_Threads.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()
  IF (Kokkos_ENABLE_OpenMP)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_OpenMP.tmpl
      "Ifpack2_CLASS_FILE_NAME_OpenMP.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_OpenMP.tmpl
      "Ifpack2_CLASS_FILE_NAME_OpenMP.cpp"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CPP_SOURCES ${IFPACK2_SRCS})
  ENDIF()
  IF (KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda)
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_CLASSES}"
      Ifpack2_ETI_SC_LO_GO_Cuda.tmpl
      "${CMAKE_CURRENT_BINARY_DIR}/Ifpack2_CLASS_FILE_NAME_Cuda.cu"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CU_SOURCES ${IFPACK2_SRCS})
    IFPACK2_PROCESS_ETI_TEMPLATE(
      "${IFPACK2_ETI_LO_GO_CLASSES}"
      Ifpack2_ETI_LO_GO_Cuda.tmpl
      "${CMAKE_CURRENT_BINARY_DIR}/Ifpack2_CLASS_FILE_NAME_Cuda.cu"
      IFPACK2_SRCS)
    LIST(APPEND IFPACK2_CU_SOURCES ${IFPACK2_SRCS})
  ENDIF()
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${IFPACK2_CPP_SOURCES} ${IFPACK2_CUDA_SOURCES}
  )

IF(KokkosClassic_ENABLE_KokkosCompat AND Kokkos_ENABLE_Cuda)
  IF( TPL_ENABLE_MPI )
    CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    ifpack2_cuda
    SOURCES ${IFPACK2_CU_SOURCES}
    CUDALIBRARY
    )
ENDIF()
