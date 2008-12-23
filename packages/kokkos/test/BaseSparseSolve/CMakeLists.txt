
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../BaseSparseMultiply)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BaseSparseSolve
  SOURCES cxx_main.cpp
          scsrmm2a5.f
  ARGS -v
  COMM serial
  )
