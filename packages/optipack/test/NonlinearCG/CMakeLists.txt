INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_BINARY_DIR}/../thyra/core/test/nonlinear/models)
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../thyra/core/test/nonlinear/models)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NonlinearCGUnitTests
  SOURCES
    NonlinearCG_UnitTests.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  DEPLIBS thyra_nonlinear_test_models
  COMM serial mpi
  ARGS
    " "   # Make sure default works
    "--nonlin-max-iters=14"   # Default solver type is FR
    "--solver-type=PR+ --nonlin-max-iters=14 --nonlin-error-tol=5e-2"
    "--solver-type=HS --nonlin-max-iters=14"
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

# NOTE: The FR-PR method does not seem to be working?
#"--solver-type=FR-PR --nonlin-max-iters=45 --nonlin-error-tol=1e-2 --nonlin-solve-tol=1e-4"
