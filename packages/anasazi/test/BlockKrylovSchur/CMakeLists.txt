INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../MVOPTester)

IF(Teuchos_ENABLE_COMPLEX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockKrylovSchur_complex_lap_test
    SOURCES cxx_main_complex_lap.cpp
    ARGS --verbose
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BlockKrylovSchur_complex_test
    SOURCES cxx_main_complex.cpp
    COMM serial mpi
    ARGS
      "--filename=mhd1280b.cua --insitu --herm"
      "--filename=mhd1280a.cua --insitu --nonherm"
      "--filename=mhd1280b.cua --exsitu --herm"
      "--filename=mhd1280a.cua --exsitu --nonherm"
    STANDARD_PASS_OUTPUT
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockKrylovSchurComplexFiles
    SOURCE_DIR ${PACKAGE_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.cua mhd1280a.cua
    EXEDEPS BlockKrylovSchur_complex_test
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)

ENDIF()
