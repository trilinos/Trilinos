INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../TriKota/Dakota/install/lib/)

# Required for Dakota
ADD_DEFINITIONS(${TRIKOTA_DAKOTA_DEFINITIONS})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTests
  SOURCES
    Piro_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    MockModelEval_C.cpp
    MockModelEval_C.hpp
    MockModelEval_D.cpp
    MockModelEval_D.hpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_MatrixFreeOperator_UnitTests
  SOURCES
    Piro_Epetra_MatrixFreeOperator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_EpetraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MatrixFreeDecorator_UnitTests
  SOURCES
    Piro_MatrixFreeDecorator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EvalModel
  SOURCES
  Main_EvalModel.cpp
  MockModelEval_A.cpp
  MockModelEval_A.hpp
  NUM_MPI_PROCS 1-4
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

IF (Piro_ENABLE_NOX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    EpetraSolver
    SOURCES
    Main_EpetraSolver.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    SaveEigenData_Epetra.cpp
    SaveEigenData_Epetra.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraSolver
    SOURCES
    Main_ThyraSolver.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AnalysisDriver
    SOURCES
    Main_AnalysisDriver.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    ObserveSolution_Epetra.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SecondOrderIntegrator
    SOURCES
    Main_SecondOrderIntegrator.cpp
    MockModelEval_B.cpp
    MockModelEval_B.hpp
    ObserveSolution_Epetra.hpp
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NOXSolver_UnitTests
    SOURCES
    Piro_NOXSolver_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    Piro_Test_ThyraSupport.hpp
    Piro_Test_WeakenedModelEvaluator.hpp
    Piro_Test_MockObserver.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
    DEST_FILES   input_Solve_NOX_1.xml
                 input_Solve_NOX_2.xml
                 input_Solve_NOX_3.xml
		 input_Solve_NOX_Adjoint.xml
                 input_Solve_LOCA_1.xml
                 input_Solve_Rythmos_1.xml
                 input_Solve_Rythmos_1new.xml
                 input_Solve_Rythmos_2.xml
                 input_Solve_VV.xml
                 input_Solve_TR.xml
          input_Analysis_Dakota.xml dak.in
          input_Analysis_MOOCHO.xml input_Analysis_OptiPack.xml
          input_Solve_RythmosMatrixFree_Fails.xml
	  input_SGSolve.xml
	  input_problem1.xml input_problem2.xml input_coupled.xml
	  input_problem1_sg.xml input_problem2_sg.xml input_coupled_sg.xml
	  input_SGAnalysis.xml dakota_sg.in dakota_sg_target.txt
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS EpetraSolver ThyraSolver AnalysisDriver UnitTests NOXSolver_UnitTests
  )

ENDIF()

IF (Piro_ENABLE_Rythmos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RythmosSolver_UnitTests
    SOURCES
    Piro_RythmosSolver_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_RythmosSolver_UnitTests
    SOURCES
    Piro_Epetra_RythmosSolver_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_EpetraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )
ENDIF (Piro_ENABLE_Rythmos)
