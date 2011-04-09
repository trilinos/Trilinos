INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

IF(Teuchos_ENABLE_COMPLEX)

  INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../MVOPTester)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOBPCG_complex_lap_test
    SOURCES cxx_main_complex_lap.cpp
    ARGS --verbose
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOBPCG_complex_test
    SOURCES cxx_main_complex.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=mhd1280b.cua"
    STANDARD_PASS_OUTPUT
    )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestLOBPCGComplexFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/anasazi/testmatrices
    SOURCE_FILES mhd1280b.cua
    EXEDEPS LOBPCG_complex_test
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)

ENDIF()
