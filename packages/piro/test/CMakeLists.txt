
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
    NOXSolve
    SOURCES 
    Main_NOXSolve.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
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
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
    DEST_FILES   input_1.xml input_2.xml dak.in
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS NOXSolve AnalysisDriver
  )

ENDIF()

