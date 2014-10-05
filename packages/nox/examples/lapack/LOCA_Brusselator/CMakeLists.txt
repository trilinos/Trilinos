
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Brusselator_Continuation_Example
    SOURCES
    BrusselatorProblemInterface.H
    BrusselatorProblemInterface.C
    BrusselatorContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All examples passed"
    )

  IF(Teuchos_ENABLE_COMPLEX)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_Brusselator_Hopf_Example
      SOURCES
      BrusselatorProblemInterface.H
      BrusselatorProblemInterface.C
      BrusselatorHopfContinuation.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "All examples passed"
      )
  ENDIF()

ENDIF()
