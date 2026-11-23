TRIBITS_ADD_EXECUTABLE(
  FDAssembly
  SOURCES finite_difference_assembly.cpp
  COMM serial mpi
  CATEGORIES BASIC
)

TRIBITS_ADD_TEST(
  FDAssembly
  NAME FDAssembly_5Point
  ARGS "--nx=18 --ny=14"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)
