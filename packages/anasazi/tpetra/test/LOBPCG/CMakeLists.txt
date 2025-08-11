
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_test
  SOURCES cxx_main.cpp
  ARGS
    "--verbose"
    "--verbose --ortho=ICGS"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_auxtest
  SOURCES cxx_main_auxtest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_LOBPCG_simpletest
  SOURCES cxx_simple.cpp
  ARGS
    "--verbose"
  COMM serial mpi
  )



