#1 RANK
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  vector_random_fs
  SOURCES cxx_main.cpp
  NAME vector_random_fs_npr1_rhs1
  NUM_MPI_PROCS 1
  ARGS "3501 1 1 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs4
  NUM_MPI_PROCS 1
  ARGS "3501 1 1 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs5
  NUM_MPI_PROCS 1
  ARGS "3501 1 1 5"
  COMM mpi
  )

#2 RANKS
TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs1
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs4
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs5
  NUM_MPI_PROCS 2
  ARGS "3501 1 2 5"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs1
  NUM_MPI_PROCS 2
  ARGS "3501 2 2 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs4
  NUM_MPI_PROCS 2
  ARGS "3501 2 2 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs5
  NUM_MPI_PROCS 2
  ARGS "3501 2 2 5"
  COMM mpi
  )

#3 RANKS
TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs1
  NUM_MPI_PROCS 3
  ARGS "3501 1 3 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs4
  NUM_MPI_PROCS 3
  ARGS "3501 1 3 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs5
  NUM_MPI_PROCS 3
  ARGS "3501 1 3 5"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr3_rhs1
  NUM_MPI_PROCS 3
  ARGS "3501 3 3 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr3_rhs4
  NUM_MPI_PROCS 3
  ARGS "3501 3 3 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr3_rhs5
  NUM_MPI_PROCS 3
  ARGS "3501 3 3 5"
  COMM mpi
  )

#4 RANKS
TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs1
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs4
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr1_rhs5
  NUM_MPI_PROCS 4
  ARGS "3501 1 4 5"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs1
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs4
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr2_rhs5
  NUM_MPI_PROCS 4
  ARGS "3501 2 4 5"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr4_rhs1
  NUM_MPI_PROCS 4
  ARGS "3501 4 4 1"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr4_rhs4
  NUM_MPI_PROCS 4
  ARGS "3501 4 4 4"
  COMM mpi
  )

TRIBITS_ADD_TEST(
  vector_random_fs
  NAME vector_random_fs_npr4_rhs5
  NUM_MPI_PROCS 4
  ARGS "3501 4 4 5"
  COMM mpi
  )
