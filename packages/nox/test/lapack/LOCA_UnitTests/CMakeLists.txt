

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../LOCA_TestProblems)

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK AND NOX_ENABLE_LOCA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_BorderedSolve
    SOURCES BorderedSolve.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    LOCA_TransposeBorderedSolve
    SOURCES TransposeBorderedSolve.C
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS noxtestutils localapacktestproblems
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF (NOX_ENABLE_Anasazi)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_AnasaziJacobianInverse
      SOURCES AnasaziJacobianInverse.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      DEPLIBS noxtestutils localapacktestproblems
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCA_AnasaziNotConverged
      SOURCES AnasaziNotConverged.C
      COMM serial mpi
      NUM_MPI_PROCS 1
      DEPLIBS noxtestutils localapacktestproblems
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

  ENDIF()

ENDIF()
