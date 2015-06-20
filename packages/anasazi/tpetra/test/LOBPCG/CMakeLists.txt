

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_LOBPCG_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_LOBPCG_ComplexCopyFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS Tpetra_LOBPCG_Complex_test
    )

  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   Tpetra_LOBPCG_DoubleFloat_test
  #   SOURCES cxx_main_doublefloat.cpp 
  #   ARGS 
  #   COMM serial mpi
  #   )

  # TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_LOBPCG_DoubleFloatCopyFiles
  #   SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
  #   SOURCE_FILES bcsstk14.hb
  #   EXEDEPS Tpetra_LOBPCG_DoubleFloat_test
  #   )

ENDIF()

