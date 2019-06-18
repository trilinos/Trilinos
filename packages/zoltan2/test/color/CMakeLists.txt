
IF (${PROJECT_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Coloring1
    SOURCES coloring1.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_distributed
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--inputFile=simple --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_distributed
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--inputFile=commanche_dual --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_distributed
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--inputFile=complete --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_distributed
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
    "--inputFile=noedges --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_distributed_one-gets-none
  NUM_MPI_PROCS 9
  COMM serial mpi
  ARGS
    "--inputFile=complete --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_color_tests
  SOURCE_FILES
    simple.mtx
    commanche_dual.mtx
    complete.mtx
    noedges.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
