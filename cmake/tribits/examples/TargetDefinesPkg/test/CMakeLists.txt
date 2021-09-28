# No defines, just default
tribits_add_executable_and_test( testcasedefault1
  SOURCES TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns DEFAULT"
  )

tribits_add_executable_and_test( testcase1
  SOURCES TargetDefinesPkg_test.cpp
  TARGET_DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_1\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_1"
  )

tribits_add_executable_and_test( testcase2
  SOURCES TargetDefinesPkg_test.cpp
  TARGET_DEFINES -DTARGETDEFINESPKG_RETURN=\"TEST_CASE_2\"
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns TEST_CASE_2"
  )

# No defines, just default
tribits_add_executable_and_test( testcasedefault2
  SOURCES TargetDefinesPkg_test.cpp
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION
    "TargetDefinesPkg_test returns DEFAULT"
  )

add_subdirectory(deprecated_defines)
