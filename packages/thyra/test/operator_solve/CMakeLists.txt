INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE(
  test_linear_op_with_solve
  SOURCES test_linear_op_with_solve.cpp
  COMM serial mpi
  )
PACKAGE_ADD_TEST(
  test_linear_op_with_solve
  NAME test_linear_op_with_solve_n1_n2
  ARGS "--n=1" "--n=2"
  NUM_MPI_PROCS 1
  COMM serial mpi
  )
PACKAGE_ADD_TEST(
  test_linear_op_with_solve
  NAME test_linear_op_with_solve_n4
  ARGS "--n=4"
  NUM_MPI_PROCS 1
  COMM serial mpi
  XHOSTTYPE Darwin 
  )
