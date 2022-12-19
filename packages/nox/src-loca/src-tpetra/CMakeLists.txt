SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-thyra)

# Base library
APPEND_SET(HEADERS
  LOCA_Tpetra_ConstraintModelEvaluator.hpp
  LOCA_Tpetra_Factory.hpp
  LOCA_BorderedSolver_TpetraHouseholder.hpp
  LOCA_Tpetra_LowRankUpdateRowMatrix.hpp
  )

APPEND_SET(SOURCES
  LOCA_Tpetra_ConstraintModelEvaluator.cpp
  LOCA_Tpetra_Factory.cpp
  LOCA_BorderedSolver_TpetraHouseholder.cpp
  LOCA_Tpetra_LowRankUpdateRowMatrix.cpp
  )

IF (NOX_ENABLE_LOCA AND NOX_ENABLE_KOKKOS_SOLVER_STACK)
  TRIBITS_ADD_LIBRARY(
    locatpetra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS locathyra loca nox
    )
ENDIF()
