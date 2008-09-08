# $Header$ 

# Required modules
#INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

# Setup include directories and library dependencies
#INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
#LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

# Test arguments for MPI
IF(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST MPITEST)
ELSE(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST)
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: This MPI logic gets duplicated in every
# CMakeLists.txt file in every package.  Instead, it should be copied into a
# single file and get reused.  NO DUPLICATION!

# Test definitions
#TRILINOS_ADD_EXECUTABLE_AND_TEST(exit ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)

#TRILINOS_ADD_EXECUTABLE_AND_TEST(BLAS ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  BLAS_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY BLAS
  )
TRILINOS_ADD_TEST(
  BLAS_test
  DIRECTORY BLAS
  NAME BLAS_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME Comm_test SOURCES Comm_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  TEUCHOS_Comm_test  
  COMM serial mpi  
  SOURCES Comm_test.cpp
  DIRECTORY Comm
  )
TRILINOS_ADD_TEST(
  TEUCHOS_Comm_test
  DIRECTORY Comm
  NAME Comm_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME SOURCES ../../test/UnitTest/Teuchos_StandardUnitTestMain.cpp DefaultMpiComm_UnitTests.cpp)
TRILINOS_ADD_EXECUTABLE(
  UnitTest_test
  COMM serial mpi  
  SOURCES Teuchos_StandardUnitTestMain.cpp ../../test/Comm/DefaultMpiComm_UnitTests.cpp
  DIRECTORY UnitTest
  )
TRILINOS_ADD_TEST(
  UnitTest_test
  DIRECTORY UnitTest
  NAME UnitTest_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(DenseMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  DenseMatrix_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY DenseMatrix
  )
TRILINOS_ADD_TEST(
  DenseMatrix_test
  DIRECTORY DenseMatrix
  NAME DenseMatrix_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(FancyOutputting NAME FancyOutputting_test ${TEST_ARGUMENTS} SOURCES FancyOutputting_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  FancyOutputting_test  
  COMM serial mpi  
  SOURCES FancyOutputting_test.cpp
  DIRECTORY FancyOutputting
  )
TRILINOS_ADD_TEST(
  FancyOutputting_test
  DIRECTORY FancyOutputting
  NAME FancyOutputting_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(Handle ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Handle_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY Handle
  )
TRILINOS_ADD_TEST(
  Handle_test
  DIRECTORY Handle
  NAME Handle_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(LAPACK ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  LAPACK_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY LAPACK
  )
TRILINOS_ADD_TEST(
  LAPACK_test
  DIRECTORY LAPACK
  NAME LAPACK_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayRCP_test ${TEST_ARGUMENTS} SOURCES ArrayRCP_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ArrayRPC_test  
  COMM serial mpi  
  SOURCES ArrayRCP_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  ArrayRPC_test
  DIRECTORY MemoryManagement
  NAME ArrayRPC_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayView_test ${TEST_ARGUMENTS} SOURCES ArrayView_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ArrayView_test  
  COMM serial mpi  
  SOURCES ArrayView_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  ArrayView_test
  DIRECTORY MemoryManagement
  NAME ArrayView-test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Array_test ${TEST_ARGUMENTS} SOURCES Array_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Array_test  
  COMM serial mpi  
  SOURCES Array_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Array_test
  DIRECTORY MemoryManagement
  NAME Array_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Ptr_test ${TEST_ARGUMENTS} SOURCES Ptr_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Ptr_test  
  COMM serial mpi  
  SOURCES Ptr_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Ptr_test
  DIRECTORY MemoryManagement
  NAME Ptr_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME RCP_test ${TEST_ARGUMENTS} SOURCES RCP_test.cpp TestClasses.cpp)
TRILINOS_ADD_EXECUTABLE(
  RCP_test  
  COMM serial mpi  
  SOURCES RCP_test.cpp TestClasses.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  RCP_test
  DIRECTORY MemoryManagement
  NAME RCP_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Tuple_test ${TEST_ARGUMENTS} SOURCES Tuple_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Tuple_test  
  COMM serial mpi  
  SOURCES Tuple_test.cpp 
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Tuple_test
  DIRECTORY MemoryManagement
  NAME Tuple_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(ScalarTraits ${TEST_ARGUMENTS} SOURCES ScalarTraits_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ScalarTraits_test  
  COMM serial mpi  
  SOURCES ScalarTraits_test.cpp 
  DIRECTORY ScalarTraits
  )
TRILINOS_ADD_TEST(
  ScalarTraits_test
  DIRECTORY ScalarTraits
  NAME ScalarTraits_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Time ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Time_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Time
  )
TRILINOS_ADD_TEST(
  Time_test
  DIRECTORY Time
  NAME Time_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Workspace ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Workspace_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Workspace
  )
TRILINOS_ADD_TEST(
  Workspace_test
  DIRECTORY Workspace
  NAME Workspace_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(XML ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  XML_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XML
  )
TRILINOS_ADD_TEST(
  XML_test
  DIRECTORY XML
  NAME XML_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(XMLParser ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  XMLParser_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XMLParser
  )
TRILINOS_ADD_TEST(
  XMLParser_test
  DIRECTORY XMLParser
  NAME XMLParser_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(dyn_cast ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  dyn_cast_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY dyn_cast
  )
TRILINOS_ADD_TEST(
  dyn_cast_test
  DIRECTORY dyn_cast
  NAME dyn_cast_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Polling MPITEST SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Polling_test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Polling
)
TRILINOS_ADD_TEST(
  Polling_test
  DIRECTORY Polling
  NAME Polling_test
  NUM_MPI_PROCS 2
  COMM mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
)

#TRILINOS_ADD_EXECUTABLE_AND_TEST(ParameterList ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  ParameterList_test 
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY ParameterList
  )
TRILINOS_ADD_TEST(
  ParameterList_test
  DIRECTORY ParameterList
  NAME ParameterList_test
  COMM serial mpi
  #HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
)
