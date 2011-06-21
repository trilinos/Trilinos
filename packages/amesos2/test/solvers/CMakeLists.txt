INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Solver_Test
  SOURCES
    Solver_Test
  ARGS "--filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/"
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

IF (${PACKAGE_NAME}_ENABLE_SuperLU)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Superlu_UnitTests
    SOURCES
      Superlu_UnitTests
      ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Superlu_Solve_Tests
    SOURCES
      Superlu_Solve_Tests
      ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_SuperLUMT)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUMTFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/amesos2/test/solvers
    SOURCE_FILES superlumt_test.xml
    EXEDEPS Solver_Test
    )
ENDIF()

PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopyFiles
  SOURCE_FILES solvers_test.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS Solver_Test
  )

