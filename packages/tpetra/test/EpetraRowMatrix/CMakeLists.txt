
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=EpetraRowMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  EpetraRowMatrix_UnitTests
  SOURCES
    EpetraRowMatrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS ${ARGS}
  STANDARD_PASS_OUTPUT
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraEpetraRowMatrixCopyFiles
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS EpetraRowMatrix_UnitTests
  )
