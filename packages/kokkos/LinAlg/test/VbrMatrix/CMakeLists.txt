INCLUDE(PackageAddExecutableAndTest)   

PACKAGE_ADD_EXECUTABLE_AND_TEST( VbrMatrix_UnitTests 
  SOURCES   
    VbrMatrix_DefaultMultiply.cpp 
    VbrMatrix_UnitTests.cpp 
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
) 
