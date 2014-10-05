
INCLUDE(TribitsAddExecutableAndTest)

IF(Amesos_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ZoltanHG
    SOURCES cxx_main.cpp
    NUM_MPI_PROCS 1-4
    COMM mpi
    )
ENDIF()
