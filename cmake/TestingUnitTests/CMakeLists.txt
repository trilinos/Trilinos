
INCLUDE(PackageAddAdvancedTest)


PACKAGE_ADD_ADVANCED_TEST( CTestFunctions_UnitTest
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS -DTrilinos_INVOKE_TESTING_UNIT_TESTS:BOOL=ON ${TRILINOS_HOME_DIR}
    PASS_REGULAR_EXPRESSION "Final UnitTests Result: PASSED"
  )

# The above test runs CMake to invoke the unit tests in
# ./UnitTests/CMakeLists.txt.  These tests run as part of the CMake
# configuration process.  The configure actually fails (currently with error
# message "Configuring incomplete, errors occurred!" but it actually passes
# because the unit tests all pass based on the above grep.
