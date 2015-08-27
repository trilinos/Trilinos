
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVectorFiller
  SOURCES MultiVectorFiller.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )


