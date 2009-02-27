
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MemoryManagement_UnitTests
  SOURCES
    RCPNodeHandle_UnitTests.cpp
    RCP_UnitTests.cpp
    Ptr_UnitTests.cpp
    Array_UnitTest_helpers.cpp
    Array_UnitTests.cpp
    ArrayView_UnitTests.cpp
    ArrayRCP_UnitTests.cpp
    ArrayConversions_UnitTests.cpp
    TestClasses.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ArrayRCP_test  
  SOURCES ArrayRCP_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ArrayView_test  
  SOURCES ArrayView_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Array_test  
  SOURCES Array_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Ptr_test  
  SOURCES Ptr_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RCP_test  
  SOURCES RCP_test.cpp TestClasses.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tuple_test  
  SOURCES Tuple_test.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  GCC_CheckedSTL_UnitTests
#  SOURCES GCC_CheckedSTL_UnitTests.cpp  
#    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS --teuchos-suppress-startup-banner
#  PASS_REGULAR_EXPRESSION
#    "error: attempt to dereference a past-the-end iterator"
#  )
