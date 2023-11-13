
TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Bug9500
      NAME Bug9500
      SOURCES Bug9500.cpp
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple_sym
      SOURCES TpetraCrsColorer.cpp
      ARGS "--file=simple --symmetric"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple_nodistrib_sym
      ARGS "--file=simple --no-distribute --symmetric"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple
      ARGS "--file=simple "
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple_nodistrib
      ARGS "--file=simple --no-distribute "
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_west0067
      ARGS "--file=west0067"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_west0067_nodistrib
      ARGS "--file=west0067 --no-distribute"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

IF (${PROJECT_NAME}_ENABLE_Galeri)

TRIBITS_ADD_TEST(
    TpetraCrsColorer
    NAME TpetraCrsColorer_galeri1
    ARGS "--xdim=100 --ydim=1 --zdim=1"
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    TpetraCrsColorer
    NAME TpetraCrsColorer_galeri2
    ARGS "--xdim=20 --ydim=15 --zdim=10"
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    TpetraCrsColorer
    NAME TpetraCrsColorer_galeri2_sym
    ARGS "--xdim=20 --ydim=15 --zdim=10 --symmetric"
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

# Test that exposes a bug in compute_transpose_graph when 
# symmetrize and symmetric parameters are false.
# I think we got lucky with the other galeri test.
TRIBITS_ADD_TEST(
    TpetraCrsColorer
    NAME TpetraCrsColorer_troublecauser
    ARGS "--xdim=10 --ydim=20 --zdim=30"
    COMM mpi
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)
ENDIF()
  
# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_adcolor_tests
  SOURCE_FILES
    simple.mtx west0067.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

