INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_serial
  SOURCES sillyCgSolve_serial.cpp
  ARGS
    "--sym-op" # Default case
    "--unsym-op --diag-scale=1.05"
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_mpi
  SOURCES sillyCgSolve_mpi.cpp
  )


IF (WIN32 AND NOT CYGWIN)
  SET(MAX_ITERS_DIM_FACTOR 20)
ELSE()
  SET(MAX_ITERS_DIM_FACTOR 10)
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_serial
  SOURCES sillyPowerMethod_serial.cpp
  ARGS --max-iters-dim-factor=${MAX_ITERS_DIM_FACTOR}
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  exampleImplicitlyComposedLinearOperators
  SOURCES exampleImplicitlyComposedLinearOperators.cpp
  ARGS
    "--verb-level=VERB_MEDIUM"
    "--verb-level=VERB_EXTREME"
  NUM_MPI_PROCS 1
  )
