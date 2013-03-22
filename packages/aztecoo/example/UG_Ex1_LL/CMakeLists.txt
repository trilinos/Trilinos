
INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UG_Ex1_LL
  SOURCES cxx_main.cpp 
  COMM serial mpi
  )
