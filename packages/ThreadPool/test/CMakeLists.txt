
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_tpi_cpp
  SOURCES test_tpi.cpp
  DEPLIBS pthread
  COMM serial
  )
