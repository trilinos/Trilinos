
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=CrsMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_UnitTests
  SOURCES
    CrsMatrix_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  CrsMatrix_InOutTests
  SOURCES
    CrsMatrix_InOutTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TpetraCrsMatrixCopyFiles
  SOURCE_FILES west0067.rua mhd1280b.cua
  EXEDEPS CrsMatrix_UnitTests
  )
