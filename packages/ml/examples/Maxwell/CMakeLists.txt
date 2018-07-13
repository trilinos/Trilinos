

TRIBITS_ADD_EXECUTABLE(
  ml_maxwell
  SOURCES ml_maxwell.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_simple_max
  SOURCES ml_simple_max.c place_holder.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ml_read_maxwell
  SOURCES ml_read_maxwell.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  ml_read_maxwell
  ARGS "Ke.mat T.mat Kn.mat"
  NUM_MPI_PROCS 4
  COMM mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ReitzingerSchoeberl_cp
  SOURCE_FILES Ke.mat Kn.mat T.mat
  )            
