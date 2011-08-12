
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../TriKota/Dakota/install/lib/)

# Required for Dakota
ADD_DEFINITIONS(${TRIKOTA_DAKOTA_DEFINITIONS})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  UnitTests
  SOURCES
    Piro_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    MockModelEval_C.cpp
    MockModelEval_C.hpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
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
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ThyraSolver
    SOURCES 
    Main_ThyraSolver.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  IF (Piro_ENABLE_Stokhos) 
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      SGAnalysisDriver
      SOURCES 
      Main_SGAnalysisDriver.cpp
      MockModelEval_C.cpp
      MockModelEval_C.hpp
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "TEST PASSED"
    )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(copySGTestInputFiles
      DEST_FILES input_SGAnalysis.xml dakota_sg.in dakota_sg_target.txt
      SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
      SOURCE_PREFIX "_"
      EXEDEPS SGAnalysisDriver
    )
  ENDIF()

  PACKAGE_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
    DEST_FILES   input_Solve_NOX_1.xml
                 input_Solve_NOX_2.xml
		 input_Solve_NOX_Adjoint.xml
                 input_Solve_LOCA_1.xml
                 input_Solve_Rythmos_1.xml
                 input_Solve_Rythmos_2.xml
                 input_Solve_VV.xml
                 input_Solve_TR.xml
          input_Analysis_Dakota.xml dak.in
          input_Analysis_MOOCHO.xml input_Analysis_OptiPack.xml 
          input_Solve_RythmosMatrixFree_Fails.xml
	  input_SGSolve.xml 
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS EpetraSolver ThyraSolver AnalysisDriver UnitTests
  )

ENDIF()

