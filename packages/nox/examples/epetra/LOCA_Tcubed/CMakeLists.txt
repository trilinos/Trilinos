
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed_Continuation_Example
    SOURCES
    Basis.H
    Basis.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    Problem_Interface.H
    Problem_Interface.C
    Example.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed_TurningPoint_Example
    SOURCES
    Basis.H
    Basis.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    Problem_Interface.H
    Problem_Interface.C
    ExampleTPContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed_PhaseTransition_Example
    SOURCES
    Basis.H
    Basis.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    Problem_Interface.H
    Problem_Interface.C
    ExamplePhaseTransContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF(TPL_ENABLE_MF)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_Tcubed_MultiContinuation_Example
      SOURCES
      Basis.H
      Basis.C
      FiniteElementProblem.H
      FiniteElementProblem.C
      Problem_Interface.H
      Problem_Interface.C
      ExampleMultiContinuation.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      LOCATcubedMultiContinuationCopyFiles
      SOURCE_FILES MFresultsDX.cfg MFresultsDX.net
      EXEDEPS LOCA_Tcubed_MultiContinuation_Example
      )

    ENDIF()

ENDIF()
