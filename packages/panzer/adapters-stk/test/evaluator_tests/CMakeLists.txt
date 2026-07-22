

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  tGatherOrientation
#  SOURCES gather_orientation.cpp ${UNIT_TEST_DRIVER}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

IF(PANZER_HAVE_EPETRA_STACK)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tGatherSolution_Epetra
    SOURCES epetra_gather_solution.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tScatterResidual_Epetra
    SOURCES epetra_scatter_residual.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tScatterDirichletResidual_Epetra
    SOURCES epetra_scatter_dirichlet_residual.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 2
    )
ENDIF(PANZER_HAVE_EPETRA_STACK)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tGatherSolution_Tpetra
  SOURCES tpetra_gather_solution.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tGatherSolution_BlockedTpetra
  SOURCES tpetra_blocked_gather_solution.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterResidual_Tpetra
  SOURCES tpetra_scatter_residual.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterResidual_BlockedTpetra
  SOURCES tpetra_blocked_scatter_residual.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterDirichletResidual_Tpetra
  SOURCES tpetra_scatter_dirichlet_residual.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterDirichletResidual_BlockedTpetra
  SOURCES tpetra_blocked_scatter_dirichlet_residual.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tBasisTimesVector
  SOURCES basis_times_vector.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tPointBasisValuesEvaluator
  SOURCES point_basis_values_evaluator.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tDomainInterface
  SOURCES domain_interface.cpp ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )
