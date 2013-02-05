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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  teuchosSubcommTest
  SOURCES teuchosSubcommTest.cpp
  ARGS
  "4"
  NUM_MPI_PROCS 4
  COMM mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  galeriBug
  SOURCES galeri_bug.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

