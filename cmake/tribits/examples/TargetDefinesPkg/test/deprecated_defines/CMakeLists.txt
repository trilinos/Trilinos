#
# Show bad/confusing behavior of deprecated DEFINES argument!
#

# No defines, just default but it is impacted by below DEFINES!
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase_deprecated_default1
  SOURCES ../TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )

# This DEFINES is overridden by the one below!
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase_deprecated_1
  SOURCES ../TargetDefinesPkg_test.cpp
  DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_1\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )

# This is the only one that matters and it affects test cases above it!
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase_deprecated_2
  SOURCES ../TargetDefinesPkg_test.cpp
  DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_2\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )

# No defines, just default but it also is affected by last DEFINES call!
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase_deprecated_default2
  SOURCES ../TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )
