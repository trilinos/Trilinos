
INCLUDE(TribitsAddExecutableAndTest)


# ##### temporarily disabling this test
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   MultiVector_UnitTests
#   SOURCES
#     MultiVector_UnitTests
#     ${TEUCHOS_STD_UNIT_TEST_MAIN}
#   COMM serial mpi
#   STANDARD_PASS_OUTPUT
# #TODO
# #  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVector_UnitTestsXpetraSpecific
  SOURCES
    MultiVector_UnitTestsXpetraSpecific
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
