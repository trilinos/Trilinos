


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple
  SOURCES simple.c
  ARGS -v
  COMM serial mpi
  XHOSTTYPE AIX
  )
