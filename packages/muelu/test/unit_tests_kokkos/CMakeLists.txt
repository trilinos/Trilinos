INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers/BlockedTransfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers/BlockedSmoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

SET(SOURCES "")

APPEND_SET(SOURCES
  Aggregates_kokkos.cpp
  AmalgamationInfo_kokkos.cpp
  CoalesceDropFactory_kokkos.cpp
  CoarseMapFactory_kokkos.cpp
  IndexManager_kokkos.cpp
  LWGraph_kokkos.cpp
  MueLu_TestHelpers_kokkos.cpp
  MueLu_UnitTests_kokkos.cpp
  NotayAggregationFactory.cpp
  SaPFactory_kokkos.cpp
  SemiCoarsenPFactory_kokkos.cpp
  StructuredAggregation_kokkos.cpp
  TentativePFactory_kokkos.cpp
  UncoupledAggregationFactory_kokkos.cpp
  UseKokkos_kokkos.cpp
  Utilities_kokkos.cpp
)

### Tests that require other Trilinos packages

# First ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Zoltan
  )


  TRIBITS_ADD_EXECUTABLE(
    UnitTests_kokkos
    SOURCES ${SOURCES}
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

    TRIBITS_ADD_TEST(
      UnitTests_kokkos
      NAME "UnitTestsTpetra_kokkos"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_kokkos
      NAME "UnitTestsTpetra_kokkos"
      ARGS "--linAlgebra=Tpetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()
