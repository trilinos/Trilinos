INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/model_factory)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

# The following is for the file write_solution_data.hpp
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/example/assembly_engine)

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  solver
  SOURCES solver.cpp ${PACKAGE_SOURCE_DIR}/adapters/stk/example/assembly_engine/write_solution_data.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(chemical_plant_design_files
  SOURCE_FILES solver_nox.xml solver_rythmos.xml clean.sh
  EXEDEPS solver
  )