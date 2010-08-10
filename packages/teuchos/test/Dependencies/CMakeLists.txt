INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	Dependencies_test
	SOURCES main.cpp
  COMM serial)
