
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    petsc_1DFiniteDifference_example
    SOURCES Example.C FiniteDifference.H FiniteDifference.C Problem_Interface.H Problem_Interface.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    ARGS -v
    PASS_REGULAR_EXPRESSION "Converged!"
    )

ENDIF()
