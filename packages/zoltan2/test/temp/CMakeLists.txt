INCLUDE(TribitsAddExecutableAndTest)

IF (${PROJECT_NAME}_ENABLE_Zoltan)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    multivectorTest
    SOURCES multivectorTest.cpp
    ARGS
    "401"
    NUM_MPI_PROCS 4
    COMM mpi serial
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  teuchosCommTest
  SOURCES teuchosCommTest.cpp
  ARGS
  "4000"
  NUM_MPI_PROCS 4
  COMM mpi
)
