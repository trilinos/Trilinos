TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

# This is an example; it builds, but does NOT run as a test.
# You may run it manually in order to see how it works.
TRIBITS_ADD_EXECUTABLE(
  SolverFactory_example
  COMM serial mpi  
  SOURCES main.cpp PackageA.cpp PackageB.cpp PackageC.cpp
  )
