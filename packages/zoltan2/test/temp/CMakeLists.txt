INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  multivectorTest
  SOURCES multivectorTest.cpp
  NUM_MPI_PROCS 4
  COMM mpi serial
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  trilinosTest
#  SOURCES trilinosTest.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi serial
#  )
