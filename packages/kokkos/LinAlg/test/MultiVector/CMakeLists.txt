INCLUDE(PackageAddExecutableAndTest)   
PACKAGE_ADD_EXECUTABLE_AND_TEST(  MultiVector_ComputeTests  
  SOURCES 
    MultiVector_ComputeTests.cpp  
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
PACKAGE_ADD_EXECUTABLE_AND_TEST( MultiVector_UnitTests 
  SOURCES   
    MultiVector_UnitTests.cpp 
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
