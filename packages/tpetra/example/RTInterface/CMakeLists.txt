INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsETISupport)

# RTI examples can't use the GPU, so any HybridPlatform-supported example 
# must be disabled when GPU is enabled

JOIN(allnodes   "|" FALSE ${Tpetra_ETI_NODES}  )

# this example uses DefaultNode (no GPU), but does use C++11
IF (${PROJECT_NAME}_ENABLE_CXX11)

  # needs float,int,int,DN and double,int,int,DN
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra "S=float  LO=int GO=int N=${KokkosClassic_DefaultNode}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra "S=double LO=int GO=int N=${KokkosClassic_DefaultNode}")

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIExample
    SOURCES RTIExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
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
    RTITutorialOperator
    SOURCES RTITutorialOperator
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

# these use HybridPlatform, so must be disabled when GPU is enabled
# they also use CXX11
IF (${PROJECT_NAME}_ENABLE_CXX11)

  TRIBITS_ETI_TYPE_EXPANSION(dfneeds "S=float|double" "LO=int" "GO=int" "N=${allnodes}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${dfneeds})

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTITutorialSimple
    SOURCES RTITutorialSimple
    ARGS 
    COMM mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RTIInlineCG
    SOURCES RTIInlineCG
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  IF (Tpetra_ENABLE_QD AND (NOT Tpetra_ENABLE_EXPLICIT_INSTANTIATION OR Tpetra_INST_QD_REAL))
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      RTIInlineCG_qd
      SOURCES RTIInlineCG_qd
      ARGS 
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )
  ENDIF()

ENDIF()
