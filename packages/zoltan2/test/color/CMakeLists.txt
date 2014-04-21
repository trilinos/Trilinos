INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Coloring1
  SOURCES coloring1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

#TRIBITS_ADD_TEST(
#  Coloring1
#  NAME Coloring1_FF
#  NUM_MPI_PROCS 4
#  COMM serial mpi
#  ARGS
#   "--color_method=firstfit"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
