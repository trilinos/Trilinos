# $Header$

# Required modules
INCLUDE(Trilinos_Add_Executable)

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
TRILINOS_ADD_EXECUTABLE(BlockMap ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)
# Problem setting up MPI version of following test due to argument requirements
TRILINOS_ADD_EXECUTABLE(BasicPerfTest TEST INSTALL SOURCES cxx_main.cpp ARGS 16 12 1 1 25 -v)
TRILINOS_ADD_EXECUTABLE(Comm ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)

#TRILINOS_ADD_EXECUTABLE(CrsGraph ${TEST_ARGUMENTS}  SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_TARGET(
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
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_TARGET(
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


#TRILINOS_ADD_EXECUTABLE(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)

TRILINOS_ADD_EXECUTABLE(RowMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(CrsRectMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(Directory ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(FECrsGraph ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE(FECrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(FEVbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE(FEVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(ImportExport ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(IntSerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(Map ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)
# This test was commented out in the AUTOTOOLS configuration:
TRILINOS_ADD_EXECUTABLE(MapColoring ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(MultiVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(Object ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(RowMatrixTransposer ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(SerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(SerialSpdDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
# Following test is not building on Windows
#TRILINOS_ADD_EXECUTABLE(VbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(Vector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
