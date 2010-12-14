
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Required for Dakota
ADD_DEFINITIONS(-DHAVE_CONFIG_H)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EvalModel
  SOURCES 
  Main_EvalModel.cpp
  MockModelEval_A.cpp
  MockModelEval_A.hpp
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

IF (Piro_ENABLE_NOX) 
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Solver
    SOURCES 
    Main_Solver.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    SaveEigenData_Epetra.cpp
    SaveEigenData_Epetra.hpp
    COMM serial mpi
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
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
    DEST_FILES   input_Solve_NOX_1.xml
                 input_Solve_LOCA_1.xml
                 input_Solve_Rythmos_1.xml
          input_Analysis_Dakota.xml dak.in
          input_Analysis_MOOCHO.xml input_Analysis_OptiPack.xml 
          input_Solve_RythmosMatrixFree_Fails.xml
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS Solver AnalysisDriver
  )

ENDIF()

