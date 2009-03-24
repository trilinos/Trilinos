INCLUDE(PackageAddTest)
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

  MESSAGE("Add a single basic test")
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



#
# Execute the unit tests
#


# Set up some global environment stuff
SET(${PROJECT_NAME}_HOSTNAME testhost.nowhere.com)


# Assume that all unit tests will pass by default
GLOBAL_SET(UNITTEST_OVERALL_PASS TRUE)
GLOBAL_SET(UNITTEST_OVERALL_NUMPASSED 0)


# Set up the PACKAGE_ADD_TEST(...) function to only capture the output and not
# actually add the test.
SET( PACKAGE_ADD_TEST_ADD_TEST_CAPTURE TRUE )
SET( PACKAGE_ADD_TEST_ADD_TEST_SKIP TRUE )


MESSAGE("\n***")
MESSAGE("*** Testing basic usage of PACAKGE_ADD_TEST(...)")
MESSAGE("***\n")

UNITEST_PACKAGE_ADD_TEST_BASIC()


MESSAGE("\n***")
MESSAGE("*** Determine final result of all unit tests")
MESSAGE("***\n")

# Pass in the number of expected tests that must pass!
UNITTEST_FINAL_RESULT(4)
