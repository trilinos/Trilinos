INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../TPL)

SET(SOURCES
  UnitTestMain.cpp 
  TestCrsMatrix.cpp
  ../../TPL/gtest/gtest-all.cc
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixUnitTest
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  DEPLIBS kokkoscore kokkoscorelinalg
  )

#  DEPLIBS kokkoscore kokkoslinalg
  
