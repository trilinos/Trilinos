
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

# Test arguments for MPI
IF(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST MPITEST)
ELSE(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST)
ENDIF(TRILINOS_ENABLE_MPI)

# Test definitions
#TRILINOS_ADD_EXECUTABLE_AND_TEST(BlockMap ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)

TRILINOS_ADD_EXECUTABLE(
  BlockMap-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp checkmap.cpp
  DIRECTORY BlockMap
  )
TRILINOS_ADD_TEST(
  BlockMap-test
  DIRECTORY BlockMap
  NAME BlockMap-test
  ARGS -v
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

# Problem setting up MPI version of following test due to argument requirements
#TRILINOS_ADD_EXECUTABLE_AND_TEST(BasicPerfTest TEST INSTALL SOURCES cxx_main.cpp ARGS 16 12 1 1 25 -v)

TRILINOS_ADD_EXECUTABLE(
  BasicPerfTest-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp
  DIRECTORY BasicPerfTest
  )
TRILINOS_ADD_TEST(
  BasicPerfTest-test
  DIRECTORY BasicPerfTest
  NAME BasicPerfTest-test
  ARGS "16 12 1 1 25 -v"
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  Comm-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp
  DIRECTORY Comm
  )
TRILINOS_ADD_TEST(
  Comm-test
  DIRECTORY Comm
  NAME Comm-test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsGraph ${TEST_ARGUMENTS}  SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  CrsGraph-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp 
  DIRECTORY CrsGraph
  )
TRILINOS_ADD_TEST(
  CrsGraph-test
  DIRECTORY CrsGraph
  KEYWORDS unit
  NAME CrsGraph-test
  ARGS -v
  COMM mpi
  HOST somehost localhost.localdomain gabriel.sandia.gov
  
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  CrsMatrix-test 
  
  COMM serial mpi
  
  SOURCES cxx_main.cpp 
  DIRECTORY CrsMatrix
  )
TRILINOS_ADD_TEST(
  CrsMatrix-test
  DIRECTORY CrsMatrix
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  NAME CrsMatrix-test-big
  ARGS -v "1 a b c d" 2 3 4 5
  PASS_REGULAR_EXPRESSION "tests PASSED"
  
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(RowMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  RowMatrix-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp
  DIRECTORY RowMatrix
  )
TRILINOS_ADD_TEST(
  RowMatrix-test
  DIRECTORY RowMatrix
  NAME RowMatrix-test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(CrsRectMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  CrsRectMatrix-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp
  DIRECTORY CrsRectMatrix
  )
TRILINOS_ADD_TEST(
  CrsRectMatrix-test
  DIRECTORY CrsRectMatrix
  NAME CrsRectMatrix-test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Directory ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  Directory-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp
  DIRECTORY Directory
  )
TRILINOS_ADD_TEST(
  Directory-test
  DIRECTORY Directory
  NAME Directory-test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )



#TRILINOS_ADD_EXECUTABLE_AND_TEST(FECrsGraph ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
TRILINOS_ADD_EXECUTABLE(
  FECrsGraph-test 
  
  COMM mpi  
 
  SOURCES cxx_main.cpp ExecuteTestProblems.cpp
  DIRECTORY FECrsGraph
  )
TRILINOS_ADD_TEST(
  FECrsGraph-test
  DIRECTORY FECrsGraph
  NAME FECrsGraph-test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
 
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )




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
