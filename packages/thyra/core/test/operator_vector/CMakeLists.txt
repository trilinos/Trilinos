INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_product_space
  SOURCES test_product_space.cpp
  )
TRIBITS_ADD_TEST(
  test_product_space
  NAME test_product_space_serial
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_scalar_product
  SOURCES test_scalar_product.cpp
  ARGS --verbose
  )
TRIBITS_ADD_TEST(
  test_scalar_product
  NAME test_scalar_product_serial
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE(
  test_composite_linear_ops
  SOURCES test_composite_linear_ops.cpp
  )
TRIBITS_ADD_TEST(
  test_composite_linear_ops
  NAME test_composite_linear_ops_serial
  ARGS "--use-serial --verbose"
  COMM serial
  )
TRIBITS_ADD_TEST(
  test_composite_linear_ops
  NAME test_composite_linear_ops_spmd
  ARGS "--use-spmd --verbose"
  COMM mpi serial
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_std_ops
  SOURCES test_std_ops.cpp
  ARGS --verbose
  )
TRIBITS_ADD_TEST(
  test_std_ops
  NAME test_std_ops_serial
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DetachedSpmdVectorViewUnitTests
  SOURCES
    DetachedSpmdVectorView_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  DetachedSpmdVectorViewUnitTests
  NAME DetachedSpmdVectorViewUnitTests_serial
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultSpmdVectorSpaceUnitTests
  SOURCES
    DefaultSpmdVectorSpace_UnitTests.cpp
    DefaultSpmdMultiVector_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultProductVectorSpaceUnitTests
  SOURCES
    DefaultProductVectorSpace_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  DefaultProductVectorSpaceUnitTests
  NAME DefaultProductVectorSpaceUnitTests_serial
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultMultiVectorProductVectorSpaceUnitTests
  SOURCES
    DefaultMultiVectorProductVectorSpace_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  DefaultMultiVectorProductVectorSpaceUnitTests
  NAME DefaultMultiVectorProductVectorSpaceUnitTests_serial
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultBlockedLinearOpUnitTests
  SOURCES
    DefaultBlockedLinearOp_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SimpleDenseLinearOpUnitTests
  SOURCES SimpleDenseLinearOp_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DefaultSpmdVectorSpace_ParallelUnitTest
  SOURCES
    DefaultSpmdVectorSpace_Parallel_UnitTest.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_TEST(
  DefaultSpmdVectorSpace_ParallelUnitTest
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 2
  )
