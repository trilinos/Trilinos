# $Header$ 

# Required modules
INCLUDE(Trilinos_Add_Executable_And_Test)

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

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
TRILINOS_ADD_EXECUTABLE_AND_TEST(BLAS ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME Comm_test SOURCES Comm_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Comm ${TEST_ARGUMENTS} NAME SOURCES ../../test/UnitTest/Teuchos_StandardUnitTestMain.cpp DefaultMpiComm_UnitTests.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(DenseMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(FancyOutputting NAME FancyOutputting_test ${TEST_ARGUMENTS} SOURCES FancyOutputting_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Handle ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(LAPACK ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayRCP_test ${TEST_ARGUMENTS} SOURCES ArrayRCP_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME ArrayView_test ${TEST_ARGUMENTS} SOURCES ArrayView_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Array_test ${TEST_ARGUMENTS} SOURCES Array_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Ptr_test ${TEST_ARGUMENTS} SOURCES Ptr_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME RCP_test ${TEST_ARGUMENTS} SOURCES RCP_test.cpp TestClasses.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(MemoryManagement NAME Tuple_test ${TEST_ARGUMENTS} SOURCES Tuple_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(ScalarTraits ${TEST_ARGUMENTS} SOURCES ScalarTraits_test.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Time ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(Workspace ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(XML ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(XMLParser ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE_AND_TEST(dyn_cast ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)

IF(TRILINOS_ENABLE_MPI)
  STRING(REGEX MATCH "[0-9]+" MPI_NPROCS ${MPI_EXECUTABLE_FLAGS})
  IF(${MPI_NPROCS} GREATER 1)
      TRILINOS_ADD_EXECUTABLE_AND_TEST(Polling MPITEST SOURCES cxx_main.cpp)
  ENDIF(${MPI_NPROCS} GREATER 1)
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: This is some pretty ugly script code to be able
# to select a range of MPI processors.  This type of logic really needs to be
# aggregated into a general macro and resused.  Also, just because we are not
# using more than once processor does not mean that we should not at least
# build the above test.

IF(${PROJECT_NAME}_ENABLE_EXTENDED)
	TRILINOS_ADD_EXECUTABLE_AND_TEST(ParameterList ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
ENDIF(${PROJECT_NAME}_ENABLE_EXTENDED)
