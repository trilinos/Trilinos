INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	Conditions_test
	SOURCES cxx_main.cpp Condition_UnitTests.cpp)
