INCLUDE(PackageAddExecutableAndTest)   

PACKAGE_ADD_EXECUTABLE_AND_TEST( VbrMatrix_UnitTests 
  SOURCES   
    VbrMatrix_DefaultOps.cpp 
    VbrMatrix_UnitTests.cpp 
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
