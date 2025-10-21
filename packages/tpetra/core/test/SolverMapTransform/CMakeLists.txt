
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SolverMapTransform_UnitTests
  SOURCES
    SolverMapTransform_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverMapTransformFileCopies
  SOURCE_FILES matA1.txt lhs1.txt rhs1.txt
  EXEDEPS SolverMapTransform_UnitTests
  )
