INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_local_HSFC_order
  SOURCES simple_local_HSFC_order.cpp
  COMM serial mpi
  )

