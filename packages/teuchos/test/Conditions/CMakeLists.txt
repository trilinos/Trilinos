INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	Conditions_test
	SOURCES cxx_main.cpp Condition_UnitTests.cpp
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
	Conditions_Serialization_test
	SOURCES cxx_main.cpp Condition_Serialization_UnitTests.cpp
  STANDARD_PASS_OUTPUT
  )
