INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultipleFactoryManager
  SOURCES MultipleFactoryManager.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )
