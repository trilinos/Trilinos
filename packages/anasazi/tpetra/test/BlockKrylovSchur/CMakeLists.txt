
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BlockKrylovSchur_ComplexLap_test
  SOURCES cxx_main_complex_lap.cpp 
  ARGS 
  COMM serial mpi
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BlockKrylovSchur_Complex_test
    SOURCES cxx_main_complex.cpp 
    ARGS 
        "--herm --verbose"
        "--nonherm --verbose"
    COMM serial mpi
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_BlockKrylovSchur_ComplexCopyFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
    SOURCE_FILES mhd1280a.cua mhd1280b.cua
    EXEDEPS Tpetra_BlockKrylovSchur_Complex_test
    )
ENDIF()

