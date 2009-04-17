
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_LOBPCG_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
    COMM serial mpi
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_LOBPCG_ComplexCopyFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS Tpetra_LOBPCG_Complex_test
    )

  # PACKAGE_ADD_EXECUTABLE_AND_TEST(
  #   Tpetra_LOBPCG_DoubleFloat_test
  #   SOURCES cxx_main_doublefloat.cpp 
  #   ARGS 
  #   COMM serial mpi
  #   )

  # PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_LOBPCG_DoubleFloatCopyFiles
  #   SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
  #   SOURCE_FILES bcsstk14.hb
  #   EXEDEPS Tpetra_LOBPCG_DoubleFloat_test
  #   )

ENDIF()

