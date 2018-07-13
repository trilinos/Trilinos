

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

IF(${PACKAGE_NAME}_ENABLE_ML)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Include)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Coarsen)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Comm)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/FEGrid)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Krylov)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/LevelWrap)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Main)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/MatrixFree)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/MLAPI)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Operator)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/RefMaxwell)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Smoother)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../ml/src/Utils)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../ml/src)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../epetra/src)
ENDIF()

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
Ifpack2_UnitTestMain.cpp
Ifpack2_UnitTestAmesos2solver.cpp
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
Ifpack2_UnitTestRILUK.cpp
Ifpack2_UnitTestTemplate.cpp
Ifpack2_UnitTestTriDiSolver.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTestsIfpack2_cp
  SOURCE_DIR ${Ifpack2_SOURCE_DIR}/test/data
  SOURCE_FILES Dfactor_bw2.mm Dfactor_bw3.mm Dfactor_bw4.mm Dfactor_bw5.mm Dfactor_bw6.mm Dfactor_bw7.mm Dfactor_bw8.mm Lfactor_bw2.mm Lfactor_bw3.mm Lfactor_bw4.mm Lfactor_bw5.mm Lfactor_bw6.mm Lfactor_bw7.mm Lfactor_bw8.mm Ufactor_bw2.mm Ufactor_bw3.mm Ufactor_bw4.mm Ufactor_bw5.mm Ufactor_bw6.mm Ufactor_bw7.mm Ufactor_bw8.mm AdditiveSchwarzMatlab_A.mm AdditiveSchwarzMatlab_rhs.mm AdditiveSchwarzMatlab_rowmap.mm AdditiveSchwarzILU_NO_sol.mm AdditiveSchwarzILU_O_sol.mm AdditiveSchwarzSGS_NO_sol.mm AdditiveSchwarzSGS_O_sol.mm
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    unit_tests
    SOURCES ${SOURCES}
    ARGS "--globally-reduce-test-result --output-show-proc-rank --output-to-root-rank-only=-1"
    COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AdditiveSchwarz
    SOURCES Ifpack2_UnitTestAdditiveSchwarz.cpp Ifpack2_UnitTestMain.cpp
    ARGS ""
    COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
		RBILUK
		SOURCES 
		Ifpack2_UnitTestRBILUK.cpp
		${TEUCHOS_STD_UNIT_TEST_MAIN}
		ARGS ""
		COMM serial mpi
		STANDARD_PASS_OUTPUT
		)

#TRIBITS_ADD_EXECUTABLE(
#    Distributed_SGS
#    SOURCES Ifpack2_DistributedRelaxationTest.cpp 
#  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ChebyshevConvergence
  SOURCES 
    Ifpack2_UnitTestChebyshev2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
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
  STANDARD_PASS_OUTPUT
  )

IF(${PACKAGE_NAME}_ENABLE_Experimental_KokkosKernels_Features)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockTriDiContainerUnitAndPerfTests
    SOURCES 
    Ifpack2_UnitTestBlockTriDiContainer.cpp
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )
ENDIF()

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_ML)
# mfh 01 Jan 2017: 
# Fix for https://github.com/trilinos/Trilinos/issues/935
ASSERT_DEFINED (${PROJECT_NAME}_ENABLE_Epetra)
IF (${PACKAGE_NAME}_ENABLE_ML AND ${PROJECT_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockRelaxationPerformance
    SOURCES
      Ifpack2_UnitTestBlockRelaxationPerf.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverFactory
  SOURCES 
  Ifpack2_SolverFactory.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalSparseTriangularSolver
  SOURCES 
    Ifpack2_UnitTestLocalSparseTriangularSolver.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalSparseTriangularSolver2
  SOURCES 
    Ifpack2_UnitTestLocalSparseTriangularSolver2.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

# This test only needs to run on 1 process.  We should already have
# separate tests for the different Container types.  (If we don't,
# that's a bug!)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ContainerFactory
  SOURCES 
    Ifpack2_UnitTestContainerFactory.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

# This test only needs to run on 1 process, since it just tests
# whether Chebyshev accepts a certain parameter of a certain type.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue234
  SOURCES 
    Ifpack2_UnitTest234.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Issue762
  SOURCES 
    762.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_Experimental_KokkosKernels_Features)
IF (${PACKAGE_NAME}_ENABLE_Experimental_KokkosKernels_Features)
  SET(SGSMT_MATRIX_FILE sherman1.mtx)

  TRIBITS_ADD_EXECUTABLE(
    SGSMT
    SOURCES 
      Ifpack2_UnitTestSGSMT.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2_UnitTest_SGSMT_CopyFile
    DEST_FILES ${SGSMT_MATRIX_FILE}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS SGSMT
  )

  TRIBITS_ADD_TEST(
    SGSMT
    NAME SGSMT_compare_with_Jacobi
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )
ENDIF ()
