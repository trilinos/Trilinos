# $Header$

# Required modules
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

# Test arguments for MPI
IF(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST MPITEST)
ELSE(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST)
ENDIF(TRILINOS_ENABLE_MPI)

# Test definitions
TRILINOS_ADD_EXECUTABLE_AND_TEST(BlockMap ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)
# Problem setting up MPI version of following test due to argument requirements
TRILINOS_ADD_EXECUTABLE_AND_TEST(BasicPerfTest TEST INSTALL SOURCES cxx_main.cpp ARGS 16 12 1 1 25 -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)

#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsGraph ${TEST_ARGUMENTS}  SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  CrsGraph_test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp 
  DIRECTORY CrsGraph
  )
TRILINOS_ADD_TEST(
  CrsGraph_test
  DIRECTORY CrsGraph
  KEYWORDS unit
  NAME CrsGraph_test
  ARGS -v
  COMM mpi
  HOST somehost localhost.localdomain
  DESEND_INTO_DIR
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  CrsMatrix_test 
  
  COMM serial mpi
  
  SOURCES cxx_main.cpp 
  DIRECTORY CrsMatrix
  )
TRILINOS_ADD_TEST(
  CrsMatrix_test
  DIRECTORY CrsMatrix
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  NAME CrsMatrix_test_big
  ARGS -v "1 a b c d" 2 3 4 5
  PASS_REGULAR_EXPRESSION "tests PASSED"
  
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)

TRILINOS_ADD_EXECUTABLE_AND_TEST(RowMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsRectMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Directory ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(FECrsGraph ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE_AND_TEST(FECrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(FEVbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE_AND_TEST(FEVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(ImportExport ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(IntSerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Map ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)
# This test was commented out in the AUTOTOOLS configuration:
TRILINOS_ADD_EXECUTABLE_AND_TEST(MapColoring ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MultiVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Object ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(RowMatrixTransposer ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(SerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(SerialSpdDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE_AND_TEST(VbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Vector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
