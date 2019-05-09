
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_IRTR_Lap_test
  SOURCES cxx_main_lap.cpp 
  ARGS 
  COMM serial mpi
  )

IF(Tpetra_INST_COMPLEX_DOUBLE)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
  IF (${PACKAGE_NAME}_ENABLE_Triutils)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_IRTR_Complex_test
      SOURCES cxx_main_complex.cpp 
      ARGS 
      COMM serial mpi
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_IRTR_ComplexCopyFiles
      SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS Tpetra_IRTR_Complex_test
      )
  ENDIF()
ENDIF()
