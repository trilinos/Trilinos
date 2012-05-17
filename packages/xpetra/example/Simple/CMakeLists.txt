INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../utils)
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Xpetra_Laplace1D-Xpetra
  SOURCES Laplace1D-Xpetra.cpp
  DEPLIBS muelu-xpetra muelu-utils muelu
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Xpetra_Laplace1D-Tpetra
  SOURCES Laplace1D-Tpetra.cpp
  DEPLIBS muelu-xpetra muelu-utils muelu
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Xpetra_Laplace1D-EZXpetra
  SOURCES Laplace1D-EZXpetra.cpp
  DEPLIBS muelu-xpetra muelu-utils muelu
  )
