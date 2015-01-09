
# As of <= 14 Mar 2012, MultiVectorFiller doesn't work yet with
# multiple MPI processes.  Thus, I'll leave it out of the MPI tests
# for now.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVectorFiller_SerialNodeTest
  SOURCES 
    MultiVectorFiller_SerialNodeTest
  COMM serial
  STANDARD_PASS_OUTPUT
  )


