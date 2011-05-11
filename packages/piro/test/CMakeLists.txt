
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/../TriKota/Dakota/install/lib/)

# Required for Dakota
ADD_DEFINITIONS(-DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H)


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

  PACKAGE_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
    DEST_FILES   input_Solve_NOX_1.xml
                 input_Solve_NOX_2.xml
                 input_Solve_LOCA_1.xml
                 input_Solve_Rythmos_1.xml
                 input_Solve_Rythmos_2.xml
                 input_Solve_VV.xml
                 input_Solve_TR.xml
          input_Analysis_Dakota.xml dak.in
          input_Analysis_MOOCHO.xml input_Analysis_OptiPack.xml 
          input_Solve_RythmosMatrixFree_Fails.xml
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS EpetraSolver ThyraSolver AnalysisDriver
  )

ENDIF()

