
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockKrylovSchur_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
    "--verbose"
    "--verbose --ortho=ICGS"
  COMM serial mpi
  )


IF (Tpetra_INST_COMPLEX_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockKrylovSchur_Complex_test
    SOURCES cxx_main_complex.cpp
    ARGS
    "--herm --verbose"
    "--nonherm --verbose"
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockKrylovSchur_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280a.cua mhd1280b.cua mhd1280a.mtx mhd1280b.mtx
    EXEDEPS Tpetra_BlockKrylovSchur_Complex_test
    )
ENDIF ()
