
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../teuchos/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/NodeAPI)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/LinAlg)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../kokkos/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../tpetra/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../tpetra/src)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
Ifpack2_UnitTestMain.cpp
Ifpack2_UnitTestTemplate.cpp
Ifpack2_UnitTestAdditiveSchwarz.cpp
Ifpack2_UnitTestBlockRelaxation.cpp
Ifpack2_UnitTestChebyshev.cpp
Ifpack2_UnitTestContainer.cpp
Ifpack2_UnitTestCreateOverlapGraph.cpp
Ifpack2_UnitTestDenseSolver.cpp
Ifpack2_UnitTestDiagonal.cpp
Ifpack2_UnitTestFactory.cpp
Ifpack2_UnitTestFiltering.cpp
Ifpack2_UnitTestHeap.cpp
Ifpack2_UnitTestIlukGraph.cpp
Ifpack2_UnitTestILUT.cpp
Ifpack2_UnitTestOverlapGraph.cpp
Ifpack2_UnitTestOverlappingRowMatrix.cpp
Ifpack2_UnitTestParameters.cpp
Ifpack2_UnitTestPartitioning.cpp
Ifpack2_UnitTestRelaxation.cpp
Ifpack2_UnitTestAmesos2solver.cpp
Ifpack2_UnitTestRILUK.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsIfpack2_cp
  SOURCE_DIR ${Ifpack2_SOURCE_DIR}/test/data
  SOURCE_FILES Dfactor_bw2.mm Dfactor_bw3.mm Dfactor_bw4.mm Dfactor_bw5.mm Dfactor_bw6.mm Dfactor_bw7.mm Dfactor_bw8.mm Lfactor_bw2.mm Lfactor_bw3.mm Lfactor_bw4.mm Lfactor_bw5.mm Lfactor_bw6.mm Lfactor_bw7.mm Lfactor_bw8.mm Ufactor_bw2.mm Ufactor_bw3.mm Ufactor_bw4.mm Ufactor_bw5.mm Ufactor_bw6.mm Ufactor_bw7.mm Ufactor_bw8.mm 
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    DEPLIBS ifpack2_tpetra_testinst
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ChebyshevConvergence
  SOURCES 
    Ifpack2_UnitTestChebyshev2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  DEPLIBS ifpack2 tpetra ifpack2_tpetra_testinst
  STANDARD_PASS_OUTPUT
  )

# This example exercises the fix for Bug 5963.  It illustrates the
# "Zero" combine mode.  This example only works with exactly 2 MPI
# processes.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AdditiveSchwarzBug5963
  SOURCES 
    Ifpack2_UnitTestAdditiveSchwarzBug5963.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM mpi
  NUM_MPI_PROCS 2
  DEPLIBS ifpack2 tpetra ifpack2_tpetra_testinst
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RILUKSingleProcessUnitTests
  SOURCES 
    Ifpack2_UnitTestSingleProcessRILUK.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  DEPLIBS ifpack2 tpetra ifpack2_tpetra_testinst
  STANDARD_PASS_OUTPUT
  )



