
INCLUDE(PackageAddExecutableAndTest)

IF (${PACKAGE_NAME}_ENABLE_COMPLEX)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    DenseMatrix_test  
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ENDIF()
