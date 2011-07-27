
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

SET(ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/")

IF (SITE STREQUAL "gabriel.sandia.gov")
  SET(ARGS
    "${ARGS} --not-unit-test=VbrMatrix_int_ComplexFloat_FullMatrixComplex_UnitTest")
ENDIF()

PACKAGE_ADD_EXECUTABLE(
  VbrMatrix_UnitTests
  SOURCES
    VbrMatrix_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetratestlib
  )

PACKAGE_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests1
  ARGS ${ARGS}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests2
  ARGS ${ARGS}
  COMM mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_TEST(
  VbrMatrix_UnitTests
  NAME VbrMatrix_UnitTests3
  ARGS ${ARGS}
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
  )

