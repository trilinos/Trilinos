INCLUDE(PackageAddExecutableAndTest)   
PACKAGE_ADD_EXECUTABLE_AND_TEST(  RelaxationTest
  SOURCES 
    RelaxationTest.cpp  
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
