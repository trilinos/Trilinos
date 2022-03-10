# Note: The tests for the Behavior class are scattered across several files,
# rather than being confined to a single file.  The reason is that the Behavior
# class is instantiated one time only and once environment variables are read,
# they are cached for future use.  Therefore, to test several values of an
# environment variable, several tests need to be created (one for each distinct
# value of the environment variable).

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Behavior_Default
  SOURCES
    Behavior_Default.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Behavior_Named
  SOURCES
    Behavior_Named.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Behavior_Off
  SOURCES
    Behavior_Off.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Behavior_On
  SOURCES
    Behavior_On.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
