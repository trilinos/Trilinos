#2 RANKS -- 2 COMMS, EACH COMM has 1 RANK
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vector_random_mc
  SOURCES cxx_main.cpp
  NAME vector_random_mc_npr1_rhs1
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr1_rhs4
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr1_rhs5
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 5"
  COMM mpi
  )

#4 RANKS -- 2 COMMS, EACH COMM has 2 RANKS
TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr1_rhs1
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr1_rhs4
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr1_rhs5
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 5"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr2_rhs1
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr2_rhs4
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_mc
  NAME vector_random_mc_npr2_rhs5
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 5"
  COMM mpi
  )
