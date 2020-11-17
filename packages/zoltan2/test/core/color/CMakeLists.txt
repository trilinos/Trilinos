
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
  NAME Coloring1_simple_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=Hybrid --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=Hybrid --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=Hybrid --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=Hybrid --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=noedges --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=complete --colorMethod=Hybrid"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=2GL --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=2GL --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=2GL --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=2GL --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=noedges --colorMethod=2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=complete --colorMethod=2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2 --prepartition=rows --verbose"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2 --prepartition=nonzeros --verbose"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=D2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D2_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=D2 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D2_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=D2 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_D2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=noedges --colorMethod=D2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_D2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=complete --colorMethod=D2"
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
