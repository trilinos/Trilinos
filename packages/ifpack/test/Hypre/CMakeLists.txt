
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  hypre_test
  SOURCES IFPACK_Hypre_UnitTest.cpp hypre_UnitTest.cpp
  COMM serial mpi
  )

