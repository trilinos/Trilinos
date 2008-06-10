# $Header$ 

INCLUDE(Trilinos_Add_Executable)

INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)

# 2008/06/09: rabartl: Above: You should not be listing include directories
# here!  The list of header files should be defined in one place and one place
# only and reused everywhere!  NO DUPLICATION!

LINK_LIBRARIES(teuchos ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})

# 2008/06/09: rabartl: Above: The link libraries required by Teuchos should
# not be defined here!  Insteady, they should be defined in one file and
# reused everwhere!  NO DUPLICATION!

IF(TRILINOS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	LINK_LIBRARIES(${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
	SET(TEST_ARGUMENTS TEST MPITEST)
ELSE(TRILINOS_ENABLE_MPI)
	SET(TEST_ARGUMENTS TEST)
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: This MPI logic gets duplicated in every
# CMakeLists.txt file in every package.  Instead, it should be copied into a
# single file and get reused.  NO DUPLICATION!

#TRILINOS_ADD_EXECUTABLE(exit ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(BLAS ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(Comm ${TEST_ARGUMENTS} SOURCES Comm_test.cpp)
TRILINOS_ADD_EXECUTABLE(DenseMatrix ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(FancyOutputting NAME FancyOutputting_test ${TEST_ARGUMENTS} SOURCES FancyOutputting_test.cpp)
TRILINOS_ADD_EXECUTABLE(Handle ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(LAPACK ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME ArrayRCP_test ${TEST_ARGUMENTS} SOURCES ArrayRCP_test.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME ArrayView_test ${TEST_ARGUMENTS} SOURCES ArrayView_test.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME Array_test ${TEST_ARGUMENTS} SOURCES Array_test.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME Ptr_test ${TEST_ARGUMENTS} SOURCES Ptr_test.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME RCP_test ${TEST_ARGUMENTS} SOURCES RCP_test.cpp)
TRILINOS_ADD_EXECUTABLE(MemoryManagement NAME Tuple_test ${TEST_ARGUMENTS} SOURCES Tuple_test.cpp)
TRILINOS_ADD_EXECUTABLE(ScalarTraits ${TEST_ARGUMENTS} SOURCES ScalarTraits_test.cpp)
TRILINOS_ADD_EXECUTABLE(Time ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(Workspace ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(XML ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(XMLParser ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(dyn_cast ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)

IF(TRILINOS_ENABLE_MPI)
  STRING(REGEX MATCH "[0-9]+" MPI_NPROCS ${MPI_EXECUTABLE_FLAGS})
  IF(${MPI_NPROCS} GREATER 1)
      TRILINOS_ADD_EXECUTABLE(Polling MPITEST SOURCES cxx_main.cpp)
  ENDIF(${MPI_NPROCS} GREATER 1)
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: This is some pretty ugly script code to be able
# to select a range of MPI processors.  This type of logic really needs to be
# aggregated into a general macro and resused.  Also, just because we are not
# using more than once processor does not mean that we should not at least
# build the above test.

IF(${PROJECT_NAME}_ENABLE_EXTENDED)
	TRILINOS_ADD_EXECUTABLE(ParameterList ${TEST_ARGUMENTS} SOURCES cxx_main.cpp)
ENDIF(${PROJECT_NAME}_ENABLE_EXTENDED)
