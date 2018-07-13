TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Allocator_UnitTest
  SOURCES
    Allocator.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

# This test exercises an atexit() hook.  This means that printing
# something in main() won't test success -- we have to test whether it
# printed something at _exit_ of main().
# 
# We have to enclose open and close parentheses characters in [],
# because CMake interprets this string as a regular expression, and
# open and close parentheses have special meaning in those
# expressions.  Furthermore, CMake's regular expressions have a
# different syntax than Perl's or POSIX's, so I can't escape them with
# a backslash.  For details, please see:
#
# http://www.cmake.org/cmake/help/v3.0/command/string.html
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Allocator_atexit
  SOURCES Allocator_atexit.cpp
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "Proc 0: Teuchos allocation tracking: Maximum allocation [(]B[)]: [0-9]*, Current allocation [(]B[)]: 0"
  )
