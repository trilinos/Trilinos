INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	Dependencies_test
	SOURCES cxx_main.cpp Dependencies_UnitTests.cpp
  STANDARD_PASS_OUTPUT
  )
