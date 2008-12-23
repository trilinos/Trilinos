
#INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)



# Test definitions
#PACKAGE_ADD_EXECUTABLE_AND_TEST(BlockMap ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)

PACKAGE_ADD_EXECUTABLE(
  BlockMap_test 
  COMM mpi serial  
  SOURCES cxx_main.cpp checkmap.cpp
  DIRECTORY BlockMap
  )
PACKAGE_ADD_TEST(
  BlockMap_test
  DIRECTORY BlockMap
  NAME BlockMap_test
  ARGS -v
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

# Problem setting up MPI version of following test due to argument requirements
#PACKAGE_ADD_EXECUTABLE_AND_TEST(BasicPerfTest TEST INSTALL SOURCES cxx_main.cpp ARGS 16 12 1 1 25 -v)

PACKAGE_ADD_EXECUTABLE(
  BasicPerfTest_test 
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY BasicPerfTest
  )
PACKAGE_ADD_TEST(
  BasicPerfTest_test
  DIRECTORY BasicPerfTest
  NAME BasicPerfTest_test
  ARGS "16 12 1 1 25 -v"
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#PACKAGE_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  Comm_test 
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY Comm
  )
PACKAGE_ADD_TEST(
  Comm_test
  DIRECTORY Comm
  NAME Comm_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#PACKAGE_ADD_EXECUTABLE_AND_TEST(CrsGraph ${TEST_ARGUMENTS}  SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  CrsGraph_test 
  COMM mpi serial  
  SOURCES cxx_main.cpp 
  DIRECTORY CrsGraph
  )
PACKAGE_ADD_TEST(
  CrsGraph_test
  DIRECTORY CrsGraph
  KEYWORDS unit
  NAME CrsGraph_test
  ARGS -v
  COMM mpi serial
  
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(CrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  CrsMatrix_test 
  COMM serial mpi
  SOURCES cxx_main.cpp 
  DIRECTORY CrsMatrix
  )
PACKAGE_ADD_TEST(
  CrsMatrix_test
  DIRECTORY CrsMatrix
  COMM serial mpi
  NUM_MPI_PROCS 1-20
  NAME CrsMatrix_test
  ARGS -v
  PASS_REGULAR_EXPRESSION "tests PASSED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(RowMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  RowMatrix_test 
  COMM mpi serial   
  SOURCES cxx_main.cpp
  DIRECTORY RowMatrix
  )
PACKAGE_ADD_TEST(
  RowMatrix_test
  DIRECTORY RowMatrix
  NAME RowMatrix_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(CrsRectMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  CrsRectMatrix_test 
  COMM mpi serial   
  SOURCES cxx_main.cpp
  DIRECTORY CrsRectMatrix
  )
PACKAGE_ADD_TEST(
  CrsRectMatrix_test
  DIRECTORY CrsRectMatrix
  NAME CrsRectMatrix_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(Directory ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  Directory_test 
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY Directory
  )
PACKAGE_ADD_TEST(
  Directory_test
  DIRECTORY Directory
  NAME Directory_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(FECrsGraph ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  FECrsGraph_test 
  COMM mpi serial   
  SOURCES cxx_main.cpp ExecuteTestProblems.cpp
  DIRECTORY FECrsGraph
  )
PACKAGE_ADD_TEST(
  FECrsGraph_test
  DIRECTORY FECrsGraph
  NAME FECrsGraph_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

# Following test is not building on Windows
#PACKAGE_ADD_EXECUTABLE_AND_TEST(FECrsMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
#PACKAGE_ADD_EXECUTABLE_AND_TEST(FEVbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  FEVbrMatrix_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp ExecuteTestProblems.cpp
  DIRECTORY FEVbrMatrix
  )
PACKAGE_ADD_TEST(
  FEVbrMatrix_test
  DIRECTORY FEVbrMatrix
  NAME FEVbrMatrix_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
# Following test is not building on Windows
#PACKAGE_ADD_EXECUTABLE_AND_TEST(FEVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(ImportExport ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  ImportExport_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY ImportExport
  )
PACKAGE_ADD_TEST(
  ImportExport_test
  DIRECTORY ImportExport
  NAME ImportExport_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(IntSerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  IntSerialDense_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY IntSerialDense
  )
PACKAGE_ADD_TEST(
  IntSerialDense_test
  DIRECTORY IntSerialDense
  NAME IntSerialDense_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(Map ${TEST_ARGUMENTS} SOURCES cxx_main.cpp checkmap.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  Map_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp checkmap.cpp
  DIRECTORY Map
  )
PACKAGE_ADD_TEST(
  Map_test
  DIRECTORY Map
  NAME Map_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

# This test was commented out in the AUTOTOOLS configuration:
#PACKAGE_ADD_EXECUTABLE_AND_TEST(MapColoring ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  MapColoring_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY MapColoring
  )
PACKAGE_ADD_TEST(
  MapColoring_test
  DIRECTORY MapColoring
  NAME MapColoring_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(MultiVector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  MultiVector_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp
  DIRECTORY MultiVector
  )
PACKAGE_ADD_TEST(
  MultiVector_test
  DIRECTORY MultiVector
  NAME MultiVector_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(Object ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  Object_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY Object
  )
PACKAGE_ADD_TEST(
  Object_test
  DIRECTORY Object
  NAME Object_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
    
#PACKAGE_ADD_EXECUTABLE_AND_TEST(RowMatrixTransposer ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  RowMatrixTransposer_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY RowMatrixTransposer
  )
PACKAGE_ADD_TEST(
  RowMatrixTransposer_test
  DIRECTORY RowMatrixTransposer
  NAME RowMatrixTransposer_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(SerialDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  SerialDense_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY SerialDense
  )
PACKAGE_ADD_TEST(
  SerialDense_test
  DIRECTORY SerialDense
  NAME SerialDense_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(SerialSpdDense ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  SerialSpdDense_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp
  DIRECTORY SerialSpdDense
  )
PACKAGE_ADD_TEST(
  SerialSpdDense_test
  DIRECTORY SerialSpdDense
  NAME SerialSpdDense_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

# Following test is not building on Windows
#PACKAGE_ADD_EXECUTABLE_AND_TEST(VbrMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp ARGS -v)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(Vector ${TEST_ARGUMENTS} SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
PACKAGE_ADD_EXECUTABLE(
  Vector_test  
  COMM mpi serial  
  SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp
  DIRECTORY Vector
  )
PACKAGE_ADD_TEST(
  Vector_test
  DIRECTORY Vector
  NAME Vector_test
  ARGS -v
  NUM_MPI_PROCS 1
  COMM mpi serial
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
