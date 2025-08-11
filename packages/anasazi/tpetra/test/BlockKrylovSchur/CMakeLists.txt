
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockKrylovSchur_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
    "--verbose"
    "--verbose --ortho=ICGS"
  COMM serial mpi
  )



IF (Tpetra_INST_COMPLEX_DOUBLE)
ENDIF ()
