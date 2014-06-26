
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../LOCA_TestProblems)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed
    SOURCES
    Tcubed.C
    COMM serial mpi
    NUM_MPI_PROCS 4
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_TcubedConstraint
    SOURCES
    Tcubed_constraint.C
    COMM serial mpi
    NUM_MPI_PROCS 4
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_TcubedTP
    SOURCES
    TcubedTP.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF(NOX_ENABLE_Stratimikos)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_TcubedTP_stratimikos
      SOURCES
      TcubedTP_stratimikos.C
      COMM serial mpi
      NUM_MPI_PROCS 2
      DEPLIBS noxtestutils locaepetratestproblems
      ARGS -v
      PASS_REGULAR_EXPRESSION "All tests passed"
      )
  ENDIF()

ENDIF()
