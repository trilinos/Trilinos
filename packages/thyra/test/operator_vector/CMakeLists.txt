INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_product_space
  SOURCES test_product_space.cpp
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_scalar_product
  SOURCES test_scalar_product.cpp
  COMM serial
  ARGS --verbose
  )


PACKAGE_ADD_EXECUTABLE(
  test_composite_linear_ops
  SOURCES test_composite_linear_ops.cpp
  )
PACKAGE_ADD_TEST(
  test_composite_linear_ops
  NAME test_composite_linear_ops_serial
  ARGS "--use-serial --verbose"
  COMM serial
  )
PACKAGE_ADD_TEST(
  test_composite_linear_ops
  NAME test_composite_linear_ops_spmd
  ARGS "--use-spmd --verbose"
  COMM mpi serial
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_std_ops
  SOURCES test_std_ops.cpp
  ARGS --verbose
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_handles
  SOURCES test_handles.cpp
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_block_op
  SOURCES test_block_op.cpp
  ARGS --verbose
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_linear_combination
  SOURCES test_linear_combination
  ARGS --verbose
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DetachedSpmdVectorViewUnitTests
  SOURCES
    DetachedSpmdVectorView_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultSpmdVectorSpaceUnitTests
  SOURCES
    DefaultSpmdVectorSpace_UnitTests.cpp
    DefaultSpmdMultiVector_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultProductVectorSpaceUnitTests
  SOURCES
    DefaultProductVectorSpace_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultBlockedLinearOpUnitTests
  SOURCES
    DefaultBlockedLinearOp_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
