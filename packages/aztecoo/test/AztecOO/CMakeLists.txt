
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RecursiveIterate
  SOURCES aztecoo_recursive_iterate.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "ALL TESTS PASSED"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RecursiveIterate2
  SOURCES aztecoo_recursive_iterate2.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "ALL TESTS PASSED"
  )

