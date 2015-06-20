

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_OrthoManagerTester
  SOURCES cxx_main.cpp 
  COMM serial mpi
  ARGS
    "--verbose --useMass --ortho=SVQB"
    "--verbose --noMass --ortho=SVQB"
    "--verbose --useMass --ortho=Basic"
    "--verbose --noMass --ortho=Basic"
    "--verbose --useMass --ortho=ICGS"
    "--verbose --noMass --ortho=ICGS"
#    "--debug"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_OrthoManagerMatTester
  SOURCES cxx_mattest.cpp 
  COMM serial mpi
  ARGS
    "--verbose --useMass --ortho=SVQB"
    "--verbose --noMass --ortho=SVQB"
    "--verbose --useMass --ortho=Basic"
    "--verbose --noMass --ortho=Basic"
    "--verbose --useMass --ortho=ICGS"
    "--verbose --noMass --ortho=ICGS"
#    "--debug"
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_OrthoManagerGenTester
  SOURCES cxx_gentest.cpp 
  COMM serial mpi
  ARGS
    "--verbose"
    "--verbose --noMass"
#    "--debug"
  STANDARD_PASS_OUTPUT
  )
