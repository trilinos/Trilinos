INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE(
  SPMD_apply_op_UnitTests
  SOURCES
    SPMD_apply_op_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  )


TRIBITS_ADD_TEST(
  SPMD_apply_op_UnitTests
  POSTFIX_AND_ARGS_0 base ""
  POSTFIX_AND_ARGS_1 dump_rtops --dump-rtops
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_TEST(
  SPMD_apply_op_UnitTests
  POSTFIX_AND_ARGS_0 parallel ""
  POSTFIX_AND_ARGS_1 parallel_dump_rtops --dump-rtops
  COMM mpi
  STANDARD_PASS_OUTPUT
  )
