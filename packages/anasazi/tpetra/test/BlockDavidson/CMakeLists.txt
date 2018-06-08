TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockDavidson_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
  COMM serial mpi
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND Tpetra_INST_COMPLEX_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockDavidson_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockDavidson_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS Tpetra_BlockDavidson_Complex_test
    )
ENDIF()
