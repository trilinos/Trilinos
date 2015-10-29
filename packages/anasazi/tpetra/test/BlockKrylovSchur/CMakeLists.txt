
IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockKrylovSchur_ComplexLap_test
    SOURCES cxx_main_complex_lap.cpp 
    ARGS 
    COMM serial mpi
    )
ENDIF()

IF (Tpetra_INST_COMPLEX_DOUBLE)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
  IF (${PACKAGE_NAME}_ENABLE_Triutils)
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
      SOURCE_FILES mhd1280a.cua mhd1280b.cua
      EXEDEPS Tpetra_BlockKrylovSchur_Complex_test
      )
  ENDIF()
ENDIF ()
