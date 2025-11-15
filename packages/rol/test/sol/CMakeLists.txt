ROL_ADD_EXECUTABLE_AND_TEST(
  solInterface
  SOURCES test_01.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  solDistributions
  SOURCES test_02.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

# moved to adapters/teuchos/tests/sol
#ROL_ADD_EXECUTABLE_AND_TEST(
#  solROLBatchManager
#  SOURCES test_03.cpp
#  ARGS PrintItAll
#  NUM_MPI_PROCS 4
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
#  ADD_DIR_TO_NAME
#  )

# moved to adapters/teuchos/tests/sol
#ROL_ADD_EXECUTABLE_AND_TEST(
#  solSROMGenerator
#  SOURCES test_04.cpp
#  ARGS PrintItAll
#  NUM_MPI_PROCS 1
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
#  ADD_DIR_TO_NAME
#  )

# moved to adapters/teuchos/tests/sol
#ROL_ADD_EXECUTABLE_AND_TEST(
#  checkParametrizedSimOptSolve
#  SOURCES test_05.cpp
#  ARGS PrintItAll
#  NUM_MPI_PROCS 1
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
#  ADD_DIR_TO_NAME
#  )

# moved to adapters/teuchos/tests/sol
#ROL_ADD_EXECUTABLE_AND_TEST(
#  checkSROMVector
#  SOURCES test_06.cpp
#  ARGS PrintItAll
#  NUM_MPI_PROCS 1
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
#  ADD_DIR_TO_NAME
#  )

# moved to adapters/teuchos/tests/sol
#ROL_ADD_EXECUTABLE_AND_TEST(
#  checkHMCR
#  SOURCES test_07.cpp
#  ARGS PrintItAll
#  NUM_MPI_PROCS 1
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
#  ADD_DIR_TO_NAME
#  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkSuperQuantileQuadrangle
  SOURCES test_08.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkSpectralRisk
  SOURCES test_09.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkAlmostSureConstraint
  SOURCES test_10.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkMeanValueConstraint
  SOURCES test_11.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkRiskNeutralConstraint
  SOURCES test_12.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkLinearRegression
  SOURCES test_13.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkProgressiveHedging
  SOURCES test_14.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkPrimalDualRisk
  SOURCES test_15.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_ADD_EXECUTABLE_AND_TEST(
  checkStochasticProblem
  SOURCES test_16.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

ROL_COPY_FILES_TO_BINARY_DIR(
  solDataCopy
  SOURCE_FILES
    input_01.xml
    input_02.xml
    input_08.xml
    input_09.xml
    input_10.xml
    input_11.xml
    input_12.xml
    input_14.xml
    input_15.xml
    input_16.xml
    points.txt
    weights.txt
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)
