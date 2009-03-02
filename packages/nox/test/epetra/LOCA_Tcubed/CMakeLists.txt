
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed
    SOURCES 
    Tcubed.C
    COMM serial mpi
    NUM_MPI_PROCS 4
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_TcubedConstraint
    SOURCES 
    Tcubed_constraint.C
    COMM serial mpi
    NUM_MPI_PROCS 4
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_TcubedTP
    SOURCES 
    TcubedTP.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
