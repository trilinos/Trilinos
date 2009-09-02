
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_BorderedSolve
    SOURCES BorderedSolve.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_TransposeBorderedSolve
    SOURCES TransposeBorderedSolve.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF (NOX_ENABLE_Anasazi)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      LOCA_AnasaziJacobianInverse
      SOURCES AnasaziJacobianInverse.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      DEPLIBS noxtestutils localapacktestproblems
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      LOCA_AnasaziNotConverged
      SOURCES AnasaziNotConverged.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      DEPLIBS noxtestutils localapacktestproblems
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

  ENDIF()
  
ENDIF()
