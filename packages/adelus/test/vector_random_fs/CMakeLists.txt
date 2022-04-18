
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vector_random_fs
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 1
  ARGS "1000 1 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NUM_MPI_PROCS 2
  ARGS "1000 1 2"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NUM_MPI_PROCS 3
  ARGS "1000 1 3"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NUM_MPI_PROCS 4
  ARGS "1000 2 4"
  COMM mpi
  )
