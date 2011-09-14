INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE(
  Verbose
  SOURCES Verbose.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )
