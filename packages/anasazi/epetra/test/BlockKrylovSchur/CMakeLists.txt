
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_solvertest
  SOURCES cxx_main_solvertest.cpp
  ARGS 
    "--verbose"
    "--verbose --ortho=ICGS"
#    "--debug"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_test
  SOURCES cxx_main.cpp
  ARGS 
    "--verbose --exsitu"
    "--verbose --insitu"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  Epetra_BKS_specialized
  SOURCES cxx_main_specialized.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_norestart_test
  SOURCES cxx_main_norestart.cpp
  ARGS 
    "--verbose"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_sym_test
  SOURCES cxx_main_sym.cpp
  ARGS 
    "--verbose"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra_BKS_nh_test
  SOURCES cxx_main_nh.cpp
  ARGS 
    "--verbose --exsitu"
    "--verbose --insitu"
#    "--debug"
#    "--debug --shortrun"
  COMM serial mpi
  )


