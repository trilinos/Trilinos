
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AztecOO_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AztecOO_bigmem_test
  SOURCES bigmem.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
