
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Blackboard
  SOURCES Blackboard.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PowerMethod
  SOURCES PowerMethod.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

IF(ML_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TwoLevelDDAdditive
    SOURCES TwoLevelDDAdditive.cpp
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )
ENDIF()

IF(ML_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiLevel
    SOURCES MultiLevel.cpp
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Richardson
  SOURCES Richardson.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )

IF(ML_ENABLE_AztecOO)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AdaptiveSA
    SOURCES AdaptiveSA.cpp
    NUM_MPI_PROCS 1-4
    COMM serial mpi
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraInterface
  SOURCES EpetraInterface.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
