

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

IF(PANZER_HAVE_EPETRA)
  SET(GATHER_SOURCES epetra_gather_solution.cpp)
  SET(SCATTER_SOURCES epetra_scatter_residual.cpp)
  SET(SCATTER_DIRICHLET_SOURCES epetra_scatter_dirichlet_residual.cpp)
ELSE(PANZER_HAVE_EPETRA)
  SET(GATHER_SOURCES tpetra_gather_solution.cpp)
  SET(SCATTER_SOURCES tpetra_scatter_residual.cpp)
  SET(SCATTER_DIRICHLET_SOURCES tpetra_scatter_dirichlet_residual.cpp)
ENDIF(PANZER_HAVE_EPETRA)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  tGatherOrientation
#  SOURCES gather_orientation.cpp ${UNIT_TEST_DRIVER}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tGatherSolution
  SOURCES ${GATHER_SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterResidual
  SOURCES ${SCATTER_SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tScatterDirichletResidual
  SOURCES ${SCATTER_DIRICHLET_SOURCES} ${UNIT_TEST_DRIVER}
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
