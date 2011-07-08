INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	RowMatrixTansposer_test
	SOURCES main.cpp)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RowMatrixTransposer_UnitTests
  SOURCES
    RowMatrixTransposer_UnitTests
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetratestlib
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TransposeFileCopies
  SOURCE_FILES a.mtx atrans.mtx
  EXEDEPS RowMatrixTransposer_UnitTests
  )

