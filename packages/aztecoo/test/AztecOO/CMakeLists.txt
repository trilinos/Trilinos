
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AztecOO_test
  SOURCES cxx_main.cpp
  ARGS -v
  COMM serial mpi
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  AztecOO_bigmem_test
#  SOURCES bigmem.cpp
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

TRIBITS_ADD_EXECUTABLE(
  RecursiveIterate
  SOURCES aztecoo_recursive_iterate.cpp
  )

TRIBITS_ADD_TEST(
  RecursiveIterate
  NUM_MPI_PROCS 1
  POSTFIX_AND_ARGS_0 gmres --solver=gmres
  POSTFIX_AND_ARGS_1 bicgstab --solver=bicgstab
  POSTFIX_AND_ARGS_2 cg --solver=cg
  PASS_REGULAR_EXPRESSION "ALL TESTS PASSED"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RecursiveIterate2
  SOURCES aztecoo_recursive_iterate2.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "ALL TESTS PASSED"
  )

