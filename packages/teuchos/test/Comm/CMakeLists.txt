
# Required modules
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  DefaultMpiComm_UnitTests
  SOURCES
    DefaultMpiComm_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  DefaultMpiComm_UnitTests  
  NAME DefaultMpiComm_UnitTests_serial
  ARGS "--teuchos-suppress-startup-banner --show-test-details=ALL --test-serial"
  COMM serial mpi
  NUM_MPI_PROCS 1
  #STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  DefaultMpiComm_UnitTests  
  NAME DefaultMpiComm_UnitTests_mpi1
  ARGS "--teuchos-suppress-startup-banner --show-test-details=ALL"
  COMM mpi
  NUM_MPI_PROCS 1
  #STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  DefaultMpiComm_UnitTests  
  NAME DefaultMpiComm_UnitTests_mpi
  ARGS "--teuchos-suppress-startup-banner --show-test-details=ALL"
  COMM mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
  )
  # 2008/07/16: rabartl: Note: We can't run on more than 4 procs or the char
  # ordinals will overflow for std::complex<double> packet type.

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Comm_test  
  SOURCES Comm_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  ReduceAllAndScatter_test  
  COMM serial mpi  
  SOURCES ReduceAllAndScatter_test.cpp
  )
