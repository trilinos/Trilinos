INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF (${PROJECT_NAME}_ENABLE_CXX11 AND NOT Tpetra_ENABLE_Thrust)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIExample
    SOURCES RTIExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS tpetrartitestlib
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIOperatorExample
    SOURCES RTIOperatorExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIInlineCG
    SOURCES RTIInlineCG
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    DEPLIBS tpetrartitestlib
  )

  IF (Tpetra_ENABLE_QD)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      RTIInlineCG_qd
      SOURCES RTIInlineCG_qd
      ARGS 
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      DEPLIBS tpetrartitestlib
    )
  ENDIF()

ENDIF()
