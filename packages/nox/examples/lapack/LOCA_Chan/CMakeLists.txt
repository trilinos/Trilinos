
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Chan_Continuation_Example
    SOURCES
    ChanProblemInterface.H
    ChanProblemInterface.C
    ChanContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All examples passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Chan_TurningPoint_Example
    SOURCES
    ChanProblemInterface.H
    ChanProblemInterface.C
    ChanTPContinuation.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All examples passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_Chan_Homotopy_Example
    SOURCES
    ChanProblemInterface.H
    ChanProblemInterface.C
    ChanHomotopy.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "All examples passed"
    )

  IF(TPL_ENABLE_MF)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_Chan_MultiContinuation_Example
      SOURCES
      ChanProblemInterface.H
      ChanProblemInterface.C
      ChanMultiContinuation.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "All examples passed"
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      LOCAChanMultiContinuationCopyFiles
      SOURCE_FILES MFresultsDX.cfg MFresultsDX.net
      EXEDEPS LOCA_Chan_MultiContinuation_Example
      )

    ENDIF()

ENDIF()
