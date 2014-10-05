INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  solver
  SOURCES solver.cpp ${UNIT_TEST_DRIVER} user_app_RythmosObserver_Epetra.hpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(solver_test
  SOURCE_FILES solver_nox.xml solver_rythmos.xml clean.sh
  EXEDEPS solver
  )
