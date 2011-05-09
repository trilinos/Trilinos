INCLUDE(PackageLibraryMacros)
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/src/Epetra)

SET(HEADERS "")
SET(SOURCES "")

APPEND_GLOB( HEADERS 
  ./src/Epetra/*.hpp
  ./src/*.hpp
)

APPEND_GLOB( SOURCES 
  ./src/Epetra/*.cpp
  ./src/*.cpp
)

PACKAGE_ADD_EXECUTABLE(
  testdriver
  SOURCES
      ${HEADERS}
      ${SOURCES}
      ./Test_Utils.cpp
      ./Test_Utils.hpp
      ./testdriver.cpp
  COMM mpi serial
  )

PACKAGE_ADD_TEST(
   testdriver
   ARGS
   NUM_MPI_PROCS 4
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

PACKAGE_ADD_TEST(
   testdriver
   ARGS
   NUM_MPI_PROCS 1
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

SET(TEST_DRIVER_DATA_FILES
  lsc_B_2.mm
  lsc_Bt_2.mm
  lsc_F_2.mm
  lsc_Qu_2.mm
  lsc_exact_2.mm
  lsc_rhs.mm
  )

IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  SET(TEST_DRIVER_DATA_FILES
    ${TEST_DRIVER_DATA_FILES}
    nsjac.mm
    nsjac_1.mm  
    )
ENDIF()

PACKAGE_COPY_FILES_TO_BINARY_DIR(testdriver_copyfiles
  SOURCE_FILES ${TEST_DRIVER_DATA_FILES}
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS testdriver
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  IterativePreconditionerFactory_test
  SOURCES
    unit_tests/tIterativePreconditionerFactory.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LU2x2InverseOp_test
  SOURCES
    unit_tests/tLU2x2InverseOp.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  StridedEpetraOperator_test
  SOURCES
    unit_tests/tStridedEpetraOperator.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  RequestInterface_test 
  SOURCES
    unit_tests/tRequestInterface.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DiagnosticLinearOp_test 
  SOURCES
    unit_tests/tDiagnosticLinearOp.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DiagonallyScaledPreconditioner 
  SOURCES
    unit_tests/tDiagonallyScaledPreconditioner.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  InverseFactoryOperator 
  SOURCES
    unit_tests/tInverseFactoryOperator.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ProbingFactory
  SOURCES
    unit_tests/tProbingFactory.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DiagonalPreconditionerFactory 
  SOURCES
    unit_tests/tDiagonalPreconditionerFactory.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  StratimikosFactory 
  SOURCES
    unit_tests/tStratimikosFactory.cpp
    ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
