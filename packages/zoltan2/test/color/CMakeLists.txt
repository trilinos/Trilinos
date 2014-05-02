INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF (${PACKAGE_NAME}_ENABLE_Experimental)
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  Coloring1
#  SOURCES coloring1.cpp
#  NUM_MPI_PROCS 4
#  COMM serial mpi
#  ARGS
#   "--inputFile=simple"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#
ENDIF()
