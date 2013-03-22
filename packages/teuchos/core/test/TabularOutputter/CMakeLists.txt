
# Required modules
INCLUDE(TribitsAddExecutableAndTest)

ASSERT_DEFINED(CMAKE_HOST_SYSTEM_NAME)
IF (NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")

  # This test does not work on MS Windows because it prints out three
  # digits for the scientific notation while the Linux GCC compiler
  # only prints out two.

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TabularOutputterUnitTests
    COMM serial mpi  
    SOURCES
      TabularOutputter_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    )

ENDIF()
