# No defines, just default
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcasedefault1
  SOURCES TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns DEFAULT"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase1
  SOURCES TargetDefinesPkg_test.cpp
  TARGET_DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_1\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_1"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST( testcase2
  SOURCES TargetDefinesPkg_test.cpp
  TARGET_DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_2\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )

# No defines, just default
TRIBITS_ADD_EXECUTABLE_AND_TEST( testcasedefault2
  SOURCES TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns DEFAULT"
  )

ADD_SUBDIRECTORY(deprecated_defines)
