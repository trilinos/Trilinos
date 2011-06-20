INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  solver
  SOURCES solver.cpp ${UNIT_TEST_DRIVER} user_app_RythmosObserver_Epetra.hpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(solver_test
  SOURCE_FILES solver_nox.xml solver_rythmos.xml clean.sh
  EXEDEPS solver
  )