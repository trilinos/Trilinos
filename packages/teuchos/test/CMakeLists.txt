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
  BLAS-test  
  COMM mpi  
  SOURCES cxx_main.cpp
  DIRECTORY BLAS
  )
TRILINOS_ADD_TEST(
  BLAS-test
  DIRECTORY BLAS
  NAME BLAS-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME Comm_test SOURCES Comm_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  TEUCHOS-Comm-test  
  COMM mpi  
  SOURCES Comm_test.cpp
  DIRECTORY Comm
  )
TRILINOS_ADD_TEST(
  TEUCHOS-Comm-test
  DIRECTORY Comm
  NAME Comm-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME SOURCES ../../test/UnitTest/Teuchos_StandardUnitTestMain.cpp DefaultMpiComm_UnitTests.cpp)
TRILINOS_ADD_EXECUTABLE(
  UnitTest 
  COMM mpi  
  SOURCES Teuchos_StandardUnitTestMain.cpp ../../test/Comm/DefaultMpiComm_UnitTests.cpp
  DIRECTORY UnitTest
  )
TRILINOS_ADD_TEST(
  UnitTest-test
  DIRECTORY UnitTest
  NAME UnitTest-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(DenseMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  DenseMatrix-test  
  COMM mpi  
  SOURCES cxx_main.cpp
  DIRECTORY DenseMatrix
  )
TRILINOS_ADD_TEST(
  DenseMatrix-test
  DIRECTORY DenseMatrix
  NAME DenseMatrix-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(FancyOutputting NAME FancyOutputting_test ${TEST_ARGUMENTS} SOURCES FancyOutputting_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  FancyOutputting-test  
  COMM mpi  
  SOURCES FancyOutputting_test.cpp
  DIRECTORY FancyOutputting
  )
TRILINOS_ADD_TEST(
  FancyOutputting-test
  DIRECTORY FancyOutputting
  NAME FancyOutputting-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )


#TRILINOS_ADD_EXECUTABLE_AND_TEST(Handle ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Handle-test  
  COMM mpi  
  SOURCES cxx_main.cpp
  DIRECTORY Handle
  )
TRILINOS_ADD_TEST(
  Handle-test
  DIRECTORY Handle
  NAME Handle-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(LAPACK ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  LAPACK-test  
  COMM mpi  
  SOURCES cxx_main.cpp
  DIRECTORY LAPACK
  )
TRILINOS_ADD_TEST(
  LAPACK-test
  DIRECTORY LAPACK
  NAME LAPACK-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayRCP_test ${TEST_ARGUMENTS} SOURCES ArrayRCP_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ArrayRPC-test  
  COMM mpi  
  SOURCES ArrayRCP_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  ArrayRPC-test
  DIRECTORY MemoryManagement
  NAME ArrayRPC-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayView_test ${TEST_ARGUMENTS} SOURCES ArrayView_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ArrayView-test  
  COMM mpi  
  SOURCES ArrayView_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  ArrayView-test
  DIRECTORY MemoryManagement
  NAME ArrayView-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Array_test ${TEST_ARGUMENTS} SOURCES Array_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Array-test  
  COMM mpi  
  SOURCES Array_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Array-test
  DIRECTORY MemoryManagement
  NAME Array-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Ptr_test ${TEST_ARGUMENTS} SOURCES Ptr_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Ptr-test  
  COMM mpi  
  SOURCES Ptr_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Ptr-test
  DIRECTORY MemoryManagement
  NAME Ptr-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME RCP_test ${TEST_ARGUMENTS} SOURCES RCP_test.cpp TestClasses.cpp)
TRILINOS_ADD_EXECUTABLE(
  RCP-test  
  COMM mpi  
  SOURCES RCP_test.cpp TestClasses.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  RCP-test
  DIRECTORY MemoryManagement
  NAME RCP-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Tuple_test ${TEST_ARGUMENTS} SOURCES Tuple_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  Tuple-test  
  COMM mpi  
  SOURCES Tuple_test.cpp 
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Tuple-test
  DIRECTORY MemoryManagement
  NAME Tuple-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(ScalarTraits ${TEST_ARGUMENTS} SOURCES ScalarTraits_test.cpp)
TRILINOS_ADD_EXECUTABLE(
  ScalarTraits-test  
  COMM mpi  
  SOURCES ScalarTraits_test.cpp 
  DIRECTORY ScalarTraits
  )
TRILINOS_ADD_TEST(
  ScalarTraits-test
  DIRECTORY ScalarTraits
  NAME ScalarTraits-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Time ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Time-test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Time
  )
TRILINOS_ADD_TEST(
  Time-test
  DIRECTORY Time
  NAME Time-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(Workspace ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  Workspace-test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Workspace
  )
TRILINOS_ADD_TEST(
  Workspace-test
  DIRECTORY Workspace
  NAME Workspace-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(XML ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  XML-test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XML
  )
TRILINOS_ADD_TEST(
  XML-test
  DIRECTORY XML
  NAME XML-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(XMLParser ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  XMLParser-test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XMLParser
  )
TRILINOS_ADD_TEST(
  XMLParser-test
  DIRECTORY XMLParser
  NAME XMLParser-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

#TRILINOS_ADD_EXECUTABLE_AND_TEST(dyn_cast ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  dyn_cast-test  
  COMM mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY dyn_cast
  )
TRILINOS_ADD_TEST(
  dyn_cast-test
  DIRECTORY dyn_cast
  NAME dyn_cast-test
  NUM_MPI_PROCS 1
  COMM mpi
  HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )

IF(TRILINOS_ENABLE_MPI)
  STRING(REGEX MATCH "[0-9]+" MPI_NPROCS ${MPI_EXECUTABLE_FLAGS})
  IF(${MPI_NPROCS} GREATER 1)
      #TRILINOS_ADD_EXECUTABLE_AND_TEST(Polling MPITEST SOURCES cxx_main.cpp)
      TRILINOS_ADD_EXECUTABLE(
        Polling-test  
        COMM mpi  
        SOURCES cxx_main.cpp 
        DIRECTORY Polling
      )
      TRILINOS_ADD_TEST(
        Polling-test
        DIRECTORY Polling
        NAME Polling-test
        NUM_MPI_PROCS 2
        COMM mpi
        HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
        FAIL_REGULAR_EXPRESSION "tests FAILED"
      )
  ENDIF(${MPI_NPROCS} GREATER 1)
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: This is some pretty ugly script code to be able
# to select a range of MPI processors.  This type of logic really needs to be
# aggregated into a general macro and resused.  Also, just because we are not
# using more than once processor does not mean that we should not at least
# build the above test.

IF(${PROJECT_NAME}_ENABLE_EXTENDED)
	#TRILINOS_ADD_EXECUTABLE_AND_TEST(ParameterList ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
	TRILINOS_ADD_EXECUTABLE(
        ParameterList-test 
        COMM mpi  
        SOURCES cxx_main.cpp 
        DIRECTORY ParameterList
      )
      TRILINOS_ADD_TEST(
        ParameterList-test
        DIRECTORY ParameterList
        NAME ParameterList-test
        COMM mpi
        HOST somehost s868953.srn.sandia.gov gabriel.sandia.gov
        FAIL_REGULAR_EXPRESSION "tests FAILED"
      )
ENDIF(${PROJECT_NAME}_ENABLE_EXTENDED)
