
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  minres
  SOURCES test_minres.cpp
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-4"
  STANDARD_PASS_OUTPUT 
  )
