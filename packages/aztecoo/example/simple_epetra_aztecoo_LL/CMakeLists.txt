
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  epetra_crsmatrix_LL
  SOURCES cxx_main.cpp
  COMM serial mpi
  )
