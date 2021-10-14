TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockGmresKokkosEx
  SOURCES BlockGmresKokkosExFile.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BlockCGKokkosEx
  SOURCES BlockCGKokkosExFile.cpp
  COMM serial mpi
  )
