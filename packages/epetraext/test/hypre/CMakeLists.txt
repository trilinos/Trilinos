

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hypre_test
  SOURCES
    EpetraExt_hypre_UnitTest.cpp hypre_UnitTest.cpp hypre_Helpers.cpp
  COMM serial mpi
  )
