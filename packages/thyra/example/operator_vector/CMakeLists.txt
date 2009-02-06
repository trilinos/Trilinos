INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_serial
  SOURCES sillyCgSolve_serial.cpp
  ARGS
    "--use-silly-cg"
    "--use-sillier-cg"
    "--unsym-op --diag-scale=1.05"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_mpi
  SOURCES sillyCgSolve_mpi.cpp
  COMM serial mpi
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_serial
  SOURCES sillyPowerMethod_serial.cpp
  COMM serial
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  exampleImplicitlyComposedLinearOperators
  SOURCES exampleImplicitlyComposedLinearOperators.cpp
  ARGS
    "--verb-level=VERB_MEDIUM"
    "--verb-level=VERB_EXTREME"
  COMM serial
  )
