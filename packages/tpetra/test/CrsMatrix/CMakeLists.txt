
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=CrsMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    CrsMatrix_UnitTests
    SOURCES
      CrsMatrix_UnitTests
      ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    ARGS ${ARGS}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
  
  PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraCrsMatrixCopyFiles
    SOURCE_FILES west0067.rua mhd1280b.cua
    EXEDEPS CrsMatrix_UnitTests
    )
ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
