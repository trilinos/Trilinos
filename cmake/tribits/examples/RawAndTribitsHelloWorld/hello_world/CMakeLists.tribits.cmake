tribits_package(HelloWorld)
tribits_add_library(hello_world_lib
  HEADERS hello_world_lib.hpp SOURCES hello_world_lib.cpp)
tribits_add_executable(hello_world NOEXEPREFIX SOURCES hello_world_main.cpp
  INSTALLABLE)
tribits_add_test(hello_world NOEXEPREFIX PASS_REGULAR_EXPRESSION "Hello World")
tribits_add_executable_and_test(unit_tests SOURCES hello_world_unit_tests.cpp
  PASS_REGULAR_EXPRESSION "All unit tests passed")
tribits_package_postprocess()
