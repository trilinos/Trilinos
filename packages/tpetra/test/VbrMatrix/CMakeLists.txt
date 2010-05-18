
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=VbrMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  VbrMatrix_UnitTests
  SOURCES
    VbrMatrix_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  ARGS ${ARGS}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetratestlib
  )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   VbrMatrix_InOutTests
#   SOURCES
#     VbrMatrix_InOutTests
#     ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
#   ARGS ${ARGS}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
#   DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
#   )

