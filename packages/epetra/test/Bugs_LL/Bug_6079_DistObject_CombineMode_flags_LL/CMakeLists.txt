

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Bug_6079_DistObject_CombineMode_flags_LL
  SOURCES cxx_main.cpp
  ARGS -v
  COMM mpi serial
  FAIL_REGULAR_EXPRESSION "tests FAILED"
  )
