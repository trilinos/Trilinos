TRIBITS_ADD_EXECUTABLE(
    Coloring1
    SOURCES coloring1.cpp
  )

IF (${PROJECT_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_TEST(
    Coloring1
    NAME Coloring1_galeri
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

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
  "--inputFile=simple --colorMethod=D1 --timing"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_recolor_degrees
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1 --recolorDegrees"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_serial_threshold_25
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1 --serialThreshold=25"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()


TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=noedges --colorMethod=D1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_D1
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=complete --colorMethod=D1"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1-2GL --timing"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_recolor_degrees
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1-2GL --recolorDegrees"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_serial_threshold_25
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1-2GL --serialThreshold=25"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1-2GL --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D1-2GL_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D1-2GL --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1-2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1-2GL --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_D1-2GL_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=commanche_dual --colorMethod=D1-2GL --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=noedges --colorMethod=D1-2GL"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_D1-2GL
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=complete --colorMethod=D1-2GL"
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
  NAME Coloring1_simple_D2_recolor_degrees
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=D2 --recolorDegrees"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2_serial_threshold_25
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2 --serialThreshold=25"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_D2_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=D2 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

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

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

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

ENDIF()

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

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_PD2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=PD2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_PD2_recolor_degrees
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
  "--inputFile=simple --colorMethod=PD2 --recolorDegrees"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_PD2_serial_threshold_25
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=PD2 --serialThreshold=25"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_PD2_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=PD2 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_simple_PD2_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple --colorMethod=PD2 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_PD2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=PD2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_PD2_rows
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=PD2 --prepartition=rows"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_commanchedual_PD2_nonzeros
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=commanche_dual --colorMethod=PD2 --prepartition=nonzeros"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_noedges_PD2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=noedges --colorMethod=PD2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Coloring1
  NAME Coloring1_complete_PD2
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=complete --colorMethod=PD2"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )


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
