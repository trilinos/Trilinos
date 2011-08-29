INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ExtendedHashtableTest
  SOURCES testExtendedHashtable.cpp
          ExtendedHashtable.hpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )
