
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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
