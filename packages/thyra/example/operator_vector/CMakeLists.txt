INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_serial
  SOURCES sillyCgSolve_serial.cpp
  ARGS
    "--verbose --use-silly-cg"
    "--verbose --use-sillier-cg"
    "--unsym-op --diag-scale=1.05 --verbose"
  COMM serial
  )


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  sillyCgSolve_mpi
  SOURCES sillyCgSolve_mpi.cpp
  COMM serial mpi
  )


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  sillyPowerMethod_serial
  SOURCES sillyPowerMethod_serial.cpp
  COMM serial
  )


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  exampleImplicitlyComposedLinearOperators
  SOURCES exampleImplicitlyComposedLinearOperators.cpp
  ARGS
    "--verb-level=VERB_MEDIUM"
    "--verb-level=VERB_EXTREME"
  COMM serial
  )
