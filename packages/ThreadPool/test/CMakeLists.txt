
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_tpi_c
  SOURCES test_main.c test_c_tpi.c test_c_dnax.c test_tpi_unit.c
  DEPLIBS pthread
  COMM serial
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  test_tpi_cpp
  SOURCES test_tpi.cpp
  DEPLIBS pthread
  COMM serial
  )
