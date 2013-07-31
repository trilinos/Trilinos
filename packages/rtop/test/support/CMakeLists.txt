INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  supportUnitTests
  SOURCES
    validate_apply_op_UnitTests.cpp
    supportUnitTestsHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  ARGS " "
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RTOpSubRangeDecoratorUnitTests
  SOURCES
    RTOpSubRangeDecorator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SparseSubVectorUnitTests
  SOURCES
    SparseSubVector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

