
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  aztec_simple
  SOURCES ml_aztec_simple.c place_holder.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  operator
  SOURCES ml_operator.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  preconditioner
  SOURCES ml_preconditioner.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  readex
  SOURCES ml_readex.c place_holder.cpp
  COMM serial
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  user_smoothing
  SOURCES ml_user_smoothing.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  XHOSTTYPE AIX
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  read_MatrixMarket
  SOURCES ml_read_MatrixMarket.cpp
  ARGS "--matrix=A.mm"
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_ml_basic_examples
  SOURCE_FILES A.mm
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
