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
  CoalesceDropFactory_kokkos.cpp
  CoarseMapFactory_kokkos.cpp
  IndexManager_kokkos.cpp
  MueLu_TestHelpers_kokkos.cpp
  MueLu_UnitTests_kokkos.cpp
  StructuredAggregation_kokkos.cpp
  #SaPFactory_kokkos.cpp
  TentativePFactory_kokkos.cpp
)

### Tests that require other Trilinos packages

# First ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Ifpack
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Zoltan
  )

IF (${PACKAGE_NAME}_ENABLE_Galeri) #TMP
  TRIBITS_ADD_EXECUTABLE(
    UnitTests_kokkos
    SOURCES ${SOURCES}
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Epetra  AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
      ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos)

    TRIBITS_ADD_TEST(
      UnitTests_kokkos
      NAME "UnitTestsEpetra_kokkos"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      UnitTests_kokkos
      NAME "UnitTestsEpetra_kokkos"
      ARGS "--linAlgebra=Epetra"
      PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()

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

ENDIF() # Galeri TMP
