TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockDavidson_test
  SOURCES cxx_main.cpp
  ARGS
    "--verbose"
    "--verbose --exsitu --locking"
    "--verbose --exsitu --nolocking"
    "--verbose --insitu --locking"
    "--verbose --insitu --nolocking"
    "--verbose --insitu --rblocks=4 --locking"
    "--verbose --exsitu --rblocks=4 --locking"
    "--verbose --exsitu --rblocks=4 --nolocking"
    "--verbose --insitu --rblocks=4 --nolocking"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockDavidson_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockDavidson_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockDavidson_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

IF (Tpetra_INST_COMPLEX_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockDavidson_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockDavidson_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua mhd1280b.mtx
    EXEDEPS Tpetra_BlockDavidson_Complex_test
    )
ENDIF()
