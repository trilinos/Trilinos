
INCLUDE(TribitsAddExecutableAndTest)

IF(${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Test_Detailed
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    )
ENDIF()
