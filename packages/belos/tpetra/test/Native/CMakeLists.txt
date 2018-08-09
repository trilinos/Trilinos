TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_CG_variants
  SOURCES CG_variants.cpp
  COMM serial mpi
  )
