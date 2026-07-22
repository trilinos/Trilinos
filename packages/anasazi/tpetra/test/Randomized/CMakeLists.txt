IF (Tpetra_INST_COMPLEX_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_Randomized_Complex_test
    SOURCES cxx_main_complex.cpp
    ARGS 
    "--verbose"
    COMM serial mpi
    )


  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_Randomized_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua mhd1280b.mtx
    EXEDEPS Tpetra_Randomized_Complex_test
    )
ENDIF()
