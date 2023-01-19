tribits_include_directories( REQUIRED_DURING_INSTALLATION_TESTING
  ${CMAKE_CURRENT_SOURCE_DIR}/inc )

tribits_add_executable_and_test(
  HelloWorldTests
  SOURCES
    SimpleCxx_HelloWorld_Tests.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

tribits_add_test(simplecxx-helloworld  NOEXEPREFIX  NOEXESUFFIX
  NAME HelloWorldProg
  DIRECTORY ${PACKAGE_BINARY_DIR}/src
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Hello World"
  )
