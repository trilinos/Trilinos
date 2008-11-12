
# Required modules
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

ADD_SUBDIRECTORY(TabularOutputter)

# Test definitions

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
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  DefaultMpiComm_UnitTests  
  DIRECTORY Comm
  COMM serial mpi  
  SOURCES
    DefaultMpiComm_UnitTests.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  )
TRILINOS_ADD_TEST(
  DefaultMpiComm_UnitTests
  DIRECTORY Comm
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  ReduceAllAndScatter_test  
  DIRECTORY Comm
  COMM serial mpi  
  SOURCES ReduceAllAndScatter_test.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  TEUCHOS_Comm_test  
  DIRECTORY Comm
  COMM serial mpi  
  SOURCES Comm_test.cpp
  )
TRILINOS_ADD_TEST(
  TEUCHOS_Comm_test
  DIRECTORY Comm
  NAME Comm_test
  NUM_MPI_PROCS 1-10
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  Containers_test  
  DIRECTORY Containers
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )
TRILINOS_ADD_TEST(
  Containers_test  
  DIRECTORY Containers
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRILINOS_ADD_EXECUTABLE(
  UnitTest_UnitTests
  DIRECTORY UnitTest
  COMM serial mpi  
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    vector_UnitTests.cpp
    Teuchos_StandardUnitTestMain.cpp
  )
TRILINOS_ADD_TEST(
  UnitTest_UnitTests
  DIRECTORY UnitTest
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  UnitTest_BadUnitTest
  DIRECTORY UnitTest
  COMM serial mpi  
  SOURCES
    Int_UnitTests.cpp
    TemplateFunc_UnitTests.cpp
    Failing_UnitTest.cpp
    Teuchos_StandardUnitTestMain.cpp
  )
TRILINOS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_final_results
  DIRECTORY UnitTest
  ARGS "--show-test-details=FINAL_RESULTS"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Int_BadAssignment_UnitTest"
  )
TRILINOS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_failed
  DIRECTORY UnitTest
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "End Result: TEST FAILED"
  )
TRILINOS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_end_result_totals
  DIRECTORY UnitTest
  ARGS "--teuchos-suppress-startup-banner"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Summary: total = 8, run = 8, passed = 6, failed = 2"
  )
TRILINOS_ADD_TEST(
  UnitTest_BadUnitTest
  NAME UnitTest_BadUnitTest_Int_BadAssignment_failed
  DIRECTORY UnitTest
  ARGS
    "--show-test-details=FINAL_RESULTS --group-name=Int --test-name=BadAssignment"
    "--show-test-details=TEST_NAMES --group-name=Int --test-name=BadAssignment"
  NUM_MPI_PROCS 1
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "[FAILED]"
  )

TRILINOS_ADD_EXECUTABLE(
  MemoryManagement_UnitTests
  DIRECTORY MemoryManagement
  COMM serial mpi  
  SOURCES
    RCPNodeHandle_UnitTests.cpp
    RCP_UnitTests.cpp
    Array_UnitTest_helpers.cpp
    Array_UnitTests.cpp
    ArrayView_UnitTests.cpp
    ArrayRCP_UnitTests.cpp
    TestClasses.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  )
TRILINOS_ADD_TEST(
  MemoryManagement_UnitTests
  DIRECTORY MemoryManagement
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

#TRILINOS_ADD_EXECUTABLE(
#  Handle_test  
#  COMM serial mpi  
#  SOURCES cxx_main.cpp
#  DIRECTORY Handle
#  )
#TRILINOS_ADD_TEST(
#  Handle_test
#  DIRECTORY Handle
#  NAME Handle_test
#  NUM_MPI_PROCS 1
#  COMM serial mpi
#  )

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
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  ArrayRCP_test  
  COMM serial mpi  
  SOURCES ArrayRCP_test.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  ArrayRCP_test
  DIRECTORY MemoryManagement
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

TRILINOS_ADD_EXECUTABLE(
  Ptr_test  
  COMM serial mpi  
  SOURCES Ptr_test.cpp TestClasses.cpp
  DIRECTORY MemoryManagement
  )
TRILINOS_ADD_TEST(
  Ptr_test
  DIRECTORY MemoryManagement
  NAME Ptr_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )
  
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
  STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

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
  )

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
  #STANDARD_PASS_OUTPUT
  )

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
  STANDARD_PASS_OUTPUT
  )

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
  )

TRILINOS_ADD_EXECUTABLE(
  Teuchos_exit_test  
  DIRECTORY exit
  COMM serial mpi  
  SOURCES cxx_main.cpp some_c_func.c 
  )
TRILINOS_ADD_TEST(
  Teuchos_exit_test  
  DIRECTORY exit
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

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
)

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
  STANDARD_PASS_OUTPUT
)

TRILINOS_ADD_EXECUTABLE(
  Polynomial_test 
  COMM serial
  SOURCES cxx_main.cpp 
  DIRECTORY Polynomial
  )
TRILINOS_ADD_TEST(
  Polynomial_test
  DIRECTORY Polynomial
  NAME Polynomial_test
  COMM serial 
  STANDARD_PASS_OUTPUT
)
