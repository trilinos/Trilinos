
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK)


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lapack_Rosenbrock_example
    SOURCES Rosenbrock.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "Converged!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    lapack_Broyden_example
    SOURCES Broyden.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "Example Passed!"
    )

  IF(NOX_BUILD_PRERELEASE)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      lapack_Linear_example
      SOURCES Linear.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "Example Passed!"
      )
  ENDIF()

ENDIF()
