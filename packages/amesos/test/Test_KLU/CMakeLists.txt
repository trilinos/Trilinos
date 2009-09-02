
INCLUDE(PackageAddExecutableAndTest)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
IF (${PACKAGE_NAME}_ENABLE_Galeri)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Test_KLU
    SOURCES cxx_main.cpp
    ARGS 
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()
