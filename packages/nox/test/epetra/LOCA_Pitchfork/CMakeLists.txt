
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Pitchfork
    SOURCES 
    Pitchfork.C
    COMM serial mpi
    NUM_MPI_PROCS 2
    DEPLIBS noxtestutils locaepetratestproblems
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

ENDIF()
