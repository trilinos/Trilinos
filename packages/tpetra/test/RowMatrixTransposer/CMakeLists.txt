INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
	RowMatrixTransposer_test
	SOURCES main.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RowMatrixTransposer_UnitTests
  SOURCES
    RowMatrixTransposer_UnitTests
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TransposeFileCopies
  SOURCE_FILES a.mtx atrans.mtx
  EXEDEPS RowMatrixTransposer_UnitTests
  )

