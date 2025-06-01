
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  block
  SOURCES block.cpp
  NUM_MPI_PROCS 4
  IMPORTEDLIBS /Users/eeprude/install/gklib/lib/libGKlib.a -lGKlib.a
  ARGS "-L/Users/eeprude/install/gklib/lib/ -lGKlib"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  kokkosBlock
  SOURCES kokkosBlock.cpp
  NUM_MPI_PROCS 4
  ARGS "-L/Users/eeprude/install/gklib/lib/ -lGKlib"
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
