#
# RCB tests done by Zoltan in nightly testing
#

INCLUDE(TribitsAddExecutableAndTest)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  zoltanCompare_2
#  SOURCES zoltanCompare.cpp
#  NUM_MPI_PROCS 2
#  COMM serial mpi
#  ARGS "--inputFile=vwgt2
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  zoltanCompare_3
  SOURCES zoltanCompare.cpp
  NUM_MPI_PROCS 3
  COMM serial mpi
  ARGS "--inputFile=bug"
  ARGS "--inputFile=drake"
  ARGS "--inputFile=onedbug"
  ARGS "--inputFile=vwgt"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
#
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  zoltanCompare_4
#  SOURCES zoltanCompare.cpp
#  NUM_MPI_PROCS 4
#  COMM serial mpi
#  ARGS "--inputFile=ewgt"
#  ARGS "--inputFile=grid20x19"
#  ARGS "--inputFile=grid20x19 --average_cuts=yes"
#  ARGS "--inputFile=grid20x19 --rectilinear_blocks=yes"
#  ARGS "--inputFile=nograph"
#  ARGS "--inputFile=simple"
#  ARGS "--inputFile=simple --average_cuts=yes"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  zoltanCompare_5
#  SOURCES zoltanCompare.cpp
#  NUM_MPI_PROCS 5
#  COMM serial mpi
#  ARGS "--inputFile=brack2_3"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  zoltanCompare_6
#  SOURCES zoltanCompare.cpp
#  NUM_MPI_PROCS 6
#  COMM serial mpi
#  ARGS "--inputFile=degenerateAA"
#  ARGS "--inputFile=degenerate"
#  ARGS "--inputFile=degenerate --rectilinear_blocks=yes"
#  ARGS "--inputFile=hammond2"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#
#
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  zoltanCompare_8
#  SOURCES zoltanCompare.cpp
#  NUM_MPI_PROCS 8
#  COMM serial mpi
#  ARGS "--inputFile=hammond"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#
