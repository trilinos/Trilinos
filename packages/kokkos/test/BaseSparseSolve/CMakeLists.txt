
INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../BaseSparseMultiply)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  BaseSparseSolve
  SOURCES cxx_main.cpp
          scsrmm2a5.f
  ARGS -v
  COMM serial
  )
