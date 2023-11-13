#INCLUDE(TribitsAddExecutableAndTest)

#TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  Xpetra_MatrixDriver
#  SOURCES MueLu_MatrixDriver.cpp
#  DEPLIBS xpetra muelu
#  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  Xpetra_MatrixDriverTests
#  SOURCES MueLu_MatrixDriverTests.cpp
#  DEPLIBS xpetra muelu
#  )

ADD_SUBDIRECTORY(Simple)
