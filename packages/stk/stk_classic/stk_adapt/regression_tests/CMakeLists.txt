
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

FILE(GLOB SOURCES *.cpp)


# An input exodus file is required to run these regression tests.

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    adapt_RegressionTest
#    SOURCES ${SOURCES}
#    COMM serial mpi
#    STANDARD_PASS_OUTPUT
#    )

