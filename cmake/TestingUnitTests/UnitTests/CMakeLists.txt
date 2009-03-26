INCLUDE(PackageAddTest)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(UnitTestHelpers)


#
# Set up unit test functions that will be called below to actually run the
# unit tests.
#
# The reason that we use functions is so that we can change varibles just
# inside of the functions that have their own variable scoping.  In that way,
# we can keep variables that are set in one unit test from affecting the
# others.
#


FUNCTION(UNITEST_PACKAGE_ADD_TEST_BASIC)

  # Needed by PACKAGE_ADD_TEST(...)
  SET(PACKAGE_NAME PackageA)

  # Used locally
  SET(EXEN SomeExec)
  SET(PACKEXEN ${PACKAGE_NAME}_${EXEN})

  MESSAGE("Add a single basic test with no arguments")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe"
    )

  MESSAGE("Add a single basic test with a single argument")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS arg1 )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe;arg1"
    )

  MESSAGE("Add a single basic test with two arguments")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "arg1 arg2" )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe;arg1;arg2"
    )

  MESSAGE("Add two tests with simple arguments")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "arg1" "arg2 arg3" )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN}_0;./${PACKEXEN}.exe;arg1;${PACKEXEN}_1;./${PACKEXEN}.exe;arg2;arg3"
    )

  MESSAGE("Add a double quoted input argument")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "--arg1=\"bob and cats\"" )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe;--arg1=\"bob and cats\""
    )

  MESSAGE("Add a double quoted with single quotes input argument")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "--arg1=\"'bob' and 'cats'\"" )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe;--arg1=\"'bob' and 'cats'\""
    )

ENDFUNCTION()


FUNCTION(UNITEST_PACKAGE_ADD_TEST_COMM)

  # Needed by PACKAGE_ADD_TEST(...)
  SET(PACKAGE_NAME PackageB)

  # Doing default serial mode

  # Used locally
  SET(EXEN SomeExec)
  SET(PACKEXEN ${PACKAGE_NAME}_${EXEN})
  SET(MPI_EXEC_MAX_NUMPROCS 5)
  SET(MPI_EXEC mpiexec)
  SET(MPI_EXEC_PRE_NUMPROCS_FLAGS --pre-num-procs-flags)
  SET(MPI_EXEC_NUMPROCS_FLAG --num-procs)
  SET(MPI_EXEC_POST_NUMPROCS_FLAGS --post-num-procs-flags)

  MESSAGE("Add a test for serial with no COMM input")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe"
    )

  MESSAGE("Add a test for MPI with no COMM input")
  SET(TPL_ENABLE_MPI ON)
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN}_MPI_${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC};${MPI_EXEC_PRE_NUMPROCS_FLAGS};${MPI_EXEC_NUMPROCS_FLAG};${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC_POST_NUMPROCS_FLAGS};./${PACKEXEN}.exe"
    )

  # Doing MPI mode
  SET(TPL_ENABLE_MPI ON)

  MESSAGE("Add a test for MPI with no COMM input but with some args")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "arg1 arg2" )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN}_MPI_${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC};${MPI_EXEC_PRE_NUMPROCS_FLAGS};${MPI_EXEC_NUMPROCS_FLAG};${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC_POST_NUMPROCS_FLAGS};./${PACKEXEN}.exe;arg1;arg2"
    )

  MESSAGE("Add a serial-only in an MPI-only build (adds no test)")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} COMM serial )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    ""
    )

  MESSAGE("Add a test for MPI with 'COMM mpi'")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} ARGS "arg1 arg2" COMM mpi )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN}_MPI_${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC};${MPI_EXEC_PRE_NUMPROCS_FLAGS};${MPI_EXEC_NUMPROCS_FLAG};${MPI_EXEC_MAX_NUMPROCS};${MPI_EXEC_POST_NUMPROCS_FLAGS};./${PACKEXEN}.exe;arg1;arg2"
    )

  # Doing serial mode
  SET(TPL_ENABLE_MPI OFF)

  MESSAGE("Add a test for serial mode with 'COMM serial'")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} COMM serial )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    "${PACKEXEN};./${PACKEXEN}.exe"
    )

  MESSAGE("Add a test for serial mode with 'COMM mpi (adds no test)")
  GLOBAL_SET(PACKAGE_ADD_TEST_ADD_TEST_INPUT)
  PACKAGE_ADD_TEST( ${EXEN} COMM mpi )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_TEST_ADD_TEST_INPUT
    ""
    )

ENDFUNCTION()


#
# PACKAGE_ADD_ADVANCED_TEST(...) UnitTest Stuff
#


SET(PackageAddAdvancedTest_basic_cmnd_1_args_0_OUTPUT
"
#
# This is a CMake script and must be run as \"cmake -P SCRIPT_NAME\"
#

#
# Variables
#

SET( TEST_CMND_0 ls )

SET(NUM_CMNDS 1)

#
# Test invocation
#

SET(CMAKE_MODULE_PATH /home/rabartl/PROJECTS/Trilinos.base/Trilinos/cmake/utils)

INCLUDE(DriveAdvancedTest)

DRIVE_ADVANCED_TEST()
"
)


FUNCTION(UNITEST_PACKAGE_ADD_ADVANCED_TEST_BASIC)

  # Needed by PACKAGE_ADD_ADVANCED_TEST(...)
  SET(PACKAGE_NAME PackageA)

  # Used locally
  SET(EXEN SomeExec)
  SET(PACKEXEN ${PACKAGE_NAME}_${EXEN})
  SET(CMNDN ls)

  MESSAGE("***\n*** Add a single basic command with no arguments\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( PackageAddAdvancedTest_basic_cmnd_1_args_0
    TEST_0 CMND ${CMNDN}
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    "${CMNDN}"
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    1
    )
  # Check that the whole scirpt was written correctly.  I don't want to check
  # this too much because it will be hard to maintain this test if the scirpt
  # generation changes in even trival ways.
  FILE(READ
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_PackageAddAdvancedTest_basic_cmnd_1_args_0.cmake"
    PACKAGE_ADD_ADVANCED_SCRIPT)
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_SCRIPT
    "${PackageAddAdvancedTest_basic_cmnd_1_args_0_OUTPUT}"
    )

  MESSAGE("***\n*** Add a single basic command with two arguments\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( PackageAddAdvancedTest_basic_cmnd_1_args_2
    TEST_0 CMND ${CMNDN} ARGS CMakeLists.txt CMakeFiles
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    "${CMNDN} \"CMakeLists.txt\" \"CMakeFiles\""
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    1
    )
  UNITTEST_FILE_REGEX(
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_PackageAddAdvancedTest_basic_cmnd_1_args_2.cmake"
    REGEX_STRINGS
      "TEST_CMND_0 ls \"CMakeLists.txt\" \"CMakeFiles\""
    )

  MESSAGE("***\n*** Add two basic commands with 1 and two arguments\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( PackageAddAdvancedTest_basic_cmnd_2_args_1_2
    TEST_0 CMND echo ARGS "Cats and Dogs"
    TEST_1 CMND ls ARGS Cats
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    2
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    "echo \"Cats and Dogs\""
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_1
    "ls \"Cats\""
    )
  UNITTEST_FILE_REGEX(
    "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_PackageAddAdvancedTest_basic_cmnd_2_args_1_2.cmake"
    REGEX_STRINGS
      "TEST_CMND_0 echo \"Cats and Dogs\""
      "TEST_CMND_1 ls \"Cats\""
      "NUM_CMNDS 2"
    )

ENDFUNCTION()


FUNCTION(UNITEST_PACKAGE_ADD_ADVANCED_TEST_COMM)

  # Needed by PACKAGE_ADD_ADVANCED_TEST(...)
  SET(PACKAGE_NAME PackageA)

  # Used locally
  SET(EXEN SomeExec)
  SET(PACKEXEN ${PACKAGE_NAME}_${EXEN})
  SET(CMNDN ls)

  SET(PACKAGE_ADD_ADVANCED_TEST_SKIP_SCRIPT TRUE)

  # Default serial mode

  MESSAGE("***\n*** Add a test with no COMM argument\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( DummyTest
    TEST_0 CMND ${CMNDN}
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    "${CMNDN}"
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    1
    )

  MESSAGE("***\n*** Add a serial test with 'COMM serial'\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( DummyTest
    TEST_0 CMND ${CMNDN}
    COMM serial
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    "${CMNDN}"
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    1
    )

  # Doing MPI mode
  SET(TPL_ENABLE_MPI ON)

  MESSAGE("***\n*** Add serial-only test with 'COMM mpi' (will not add the test)\n***")
  PACKAGE_ADD_ADVANCED_TEST_UNITTEST_RESET()
  PACKAGE_ADD_ADVANCED_TEST( DummyTest
    TEST_0 CMND ${CMNDN}
    COMM serial
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_CMND_ARRAY_0
    ""
    )
  UNITTEST_COMPARE_CONST(
    PACKAGE_ADD_ADVANCED_TEST_NUM_CMNDS
    ""
    )

ENDFUNCTION()


#
# Execute the unit tests
#


# Set up some global environment stuff
SET(${PROJECT_NAME}_HOSTNAME testhost.nowhere.com)


# Assume that all unit tests will pass by default
GLOBAL_SET(UNITTEST_OVERALL_PASS TRUE)
GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED 0)
GLOBAL_SET(UNITTEST_OVERALL_NUMRUN 0)


# Set up the PACKAGE_ADD_TEST(...) function to only capture the output and not
# actually add the test.
SET( PACKAGE_ADD_TEST_ADD_TEST_CAPTURE TRUE )
SET( PACKAGE_ADD_TEST_ADD_TEST_SKIP TRUE )


MESSAGE("\n***")
MESSAGE("*** Testing basic usage of PACAKGE_ADD_TEST(...)")
MESSAGE("***\n")

UNITEST_PACKAGE_ADD_TEST_BASIC()

UNITEST_PACKAGE_ADD_TEST_COMM()

MESSAGE("\n***")
MESSAGE("*** Testing basic usage of PACAKGE_ADD_ADVANCED_TEST(...)")
MESSAGE("***\n")

UNITEST_PACKAGE_ADD_ADVANCED_TEST_BASIC()

UNITEST_PACKAGE_ADD_ADVANCED_TEST_COMM()

MESSAGE("\n***")
MESSAGE("*** Determine final result of all unit tests")
MESSAGE("***\n")

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(31)
