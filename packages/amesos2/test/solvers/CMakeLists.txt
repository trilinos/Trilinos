INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Solver_Test
  SOURCES Solver_Test
  )

##### SuperLu Tests #####

IF (${PACKAGE_NAME}_ENABLE_SuperLU)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superlu_test.xml
    EXEDEPS Solver_Test
    )

  PACKAGE_ADD_TEST(
    Solver_Test
    NAME SuperLU_Solver_Test
    ARGS "--xml-params=superlu_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Superlu_UnitTests
    SOURCES
      Superlu_UnitTests
      ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    CATEGORIES NIGHTLY
    STANDARD_PASS_OUTPUT
    )
ENDIF()


##### SuperLU_MT Tests #####

IF (${PACKAGE_NAME}_ENABLE_SuperLUMT)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUMTFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superlumt_test.xml
    EXEDEPS Solver_Test
    )

  PACKAGE_ADD_TEST(
    Solver_Test
    NAME SuperLU_MT_Solver_Test
    ARGS "--xml-params=superlumt_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()


##### SuperLU_DIST Tests #####

IF (${PACKAGE_NAME}_ENABLE_SuperLUDist)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopySuperLUDistFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES superludist_test.xml
    EXEDEPS Solver_Test
    )

  PACKAGE_ADD_TEST(
    Solver_Test
    NAME SuperLU_DIST_Solver_Test
    ARGS "--xml-params=superludist_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM mpi
    )
ENDIF()


##### PARDISO_MKL Tests #####

IF (${PACKAGE_NAME}_ENABLE_PARDISO_MKL)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopyPardisoMKLFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES pardiso_mkl_test.xml
    EXEDEPS Solver_Test
    )

  PACKAGE_ADD_TEST(
    Solver_Test
    NAME Pardiso_MKL_Solver_Test
    ARGS "--xml-params=pardiso_mkl_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()


##### LAPACK Tests #####

IF (${PACKAGE_NAME}_ENABLE_LAPACK)
  PACKAGE_COPY_FILES_TO_BINARY_DIR(SolverTestCopyLapackFiles
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES lapack_test.xml
    EXEDEPS Solver_Test
    )

  PACKAGE_ADD_TEST(
    Solver_Test
    NAME LAPACK_Solver_Test
    ARGS "--xml-params=lapack_test.xml --filedir=${CMAKE_CURRENT_BINARY_DIR}/../matrices/ --multiple-solves --refactor"
    STANDARD_PASS_OUTPUT
    COMM serial mpi
    )
ENDIF()
