INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Epetra64_Laplacian_Driver
  SOURCES cxx_main.cpp 
  ARGS
    "--Epetra64 --method=LOBPCG --gen=25"
    "--no-Epetra64 --method=LOBPCG --gen=25"
    "--Epetra64 --method=BKS --gen=25"
    "--no-Epetra64 --method=BKS --gen=25"
  COMM serial mpi
  )

