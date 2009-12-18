
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Main_EvalModel
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
    Main_NOXSolve
    SOURCES 
    Main_NOXSolve.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    COMM serial mpi
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
 
  PACKAGE_COPY_FILES_TO_BINARY_DIR(Piro_copyTestXMLFiles
    DEST_FILES   input_1.xml
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS Main_NOXSolve
  )

ENDIF()
