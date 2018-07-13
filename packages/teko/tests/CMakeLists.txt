
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/src/Epetra)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/src/Tpetra)

SET(HEADERS "")
SET(SOURCES "")

APPEND_GLOB( HEADERS 
  ./src/Epetra/*.hpp
  ./src/Tpetra/*.hpp
  ./src/*.hpp
)

APPEND_GLOB( SOURCES 
  ./src/Epetra/*.cpp
  ./src/Tpetra/*.cpp
  ./src/*.cpp
)

SET(UNIT_TEST_DRIVER
  "${PACKAGE_SOURCE_DIR}/tests/Teko_UnitTestMain.cpp")

# For serial tests, use 1 MPI proc, one thread when threads are available
# For parallel tests, use 4 MPI procs when OpenMP not enabled; use 2 MPI procs, 2 threads when OpenMP enabled
IF(DEFINED KOKKOS_ENABLE_OPENMP)
  IF(KOKKOS_ENABLE_OPENMP)
    SET(SERIAL_ARGS_STRING "--kokkos-threads=1")
    SET(PARALLEL_ARGS_STRING "--kokkos-threads=2")
    SET(CUSTOM_MPI_PROCS "2")
  ELSE()
    SET(SERIAL_ARGS_STRING "")
    SET(PARALLEL_ARGS_STRING "")
    SET(CUSTOM_MPI_PROCS "4")
  ENDIF()
ELSE()
  SET(SERIAL_ARGS_STRING "")
  SET(PARALLEL_ARGS_STRING "")
  SET(CUSTOM_MPI_PROCS "4")
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  testdriver
  SOURCES
      ${HEADERS}
      ${SOURCES}
      ./Test_Utils.cpp
      ./Test_Utils.hpp
      ./testdriver.cpp
  COMM mpi serial
  )

TRIBITS_ADD_TEST(
   testdriver
   ARGS ${PARALLEL_ARGS_STRING}
   NUM_MPI_PROCS ${CUSTOM_MPI_PROCS}
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

TRIBITS_ADD_TEST(
   testdriver
   ARGS ${SERIAL_ARGS_STRING}
   NUM_MPI_PROCS 1
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

IF (${PACKAGE_NAME}_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D)
  SET(testdriver_TARGET_DEFINES -DTEKO_DISABLE_LSCSTABALIZED_TPETRA_ALPAH_INV_D)
ENDIF()

TRIBITS_ADD_EXECUTABLE(
  testdriver_tpetra
  SOURCES
      ${HEADERS}
      ${SOURCES}
      ./Test_Utils.cpp
      ./Test_Utils.hpp
      ./testdriver_tpetra.cpp
  TARGET_DEFINES ${testdriver_TARGET_DEFINES}
  COMM mpi serial
  )

TRIBITS_ADD_TEST(
   testdriver_tpetra
   ARGS ${PARALLEL_ARGS_STRING}
   NUM_MPI_PROCS ${CUSTOM_MPI_PROCS}
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

TRIBITS_ADD_TEST(
   testdriver_tpetra
   ARGS ${SERIAL_ARGS_STRING}
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
  tOpMat.mm
  tOpMp.mm
  tOpRhs.mm
  szdMat.mm
  szdMp.mm
  szdRHS.mm
  )

IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  SET(TEST_DRIVER_DATA_FILES
    ${TEST_DRIVER_DATA_FILES}
    nsjac.mm
    nsjac_1.mm  
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(testdriver_copyfiles
  SOURCE_FILES ${TEST_DRIVER_DATA_FILES}
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS testdriver testdriver_tpetra
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IterativePreconditionerFactory_test
  SOURCES
    unit_tests/tIterativePreconditionerFactory.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LU2x2InverseOp_test
  SOURCES
    unit_tests/tLU2x2InverseOp.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StridedEpetraOperator_test
  SOURCES
    unit_tests/tStridedEpetraOperator.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${PARALLEL_ARGS_STRING}
  NUM_MPI_PROCS ${CUSTOM_MPI_PROCS}
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RequestInterface_test 
  SOURCES
    unit_tests/tRequestInterface.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DiagnosticLinearOp_test 
  SOURCES
    unit_tests/tDiagnosticLinearOp.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DiagonallyScaledPreconditioner 
  SOURCES
    unit_tests/tDiagonallyScaledPreconditioner.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InverseFactoryOperator 
  SOURCES
    unit_tests/tInverseFactoryOperator.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ProbingFactory
  SOURCES
    unit_tests/tProbingFactory.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DiagonalPreconditionerFactory 
  SOURCES
    unit_tests/tDiagonalPreconditionerFactory.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StratimikosFactory 
  SOURCES
    unit_tests/tStratimikosFactory.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExplicitOps 
  SOURCES
    unit_tests/tExplicitOps.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ALOperator 
  SOURCES
    unit_tests/tALOperator.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ModALPreconditioner
  SOURCES
    unit_tests/tModALPreconditioner.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReorderBlocking
  SOURCES
    unit_tests/tReorderBlocking.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UtilitiesTests 
  SOURCES
    unit_tests/tUtilitiesTests.cpp
    ${UNIT_TEST_DRIVER}
  COMM serial mpi
  ARGS ${SERIAL_ARGS_STRING}
  NUM_MPI_PROCS 1
  )
