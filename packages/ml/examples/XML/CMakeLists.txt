
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ml_XML
  SOURCES ml_XML.cpp
  NUM_MPI_PROCS 1-4
  COMM serial mpi
  )
