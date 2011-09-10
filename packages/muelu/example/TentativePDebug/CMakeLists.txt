INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#######################################################################################
#######################################################################################

PACKAGE_ADD_EXECUTABLE(
  tentativePtest
  SOURCES tentativePtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )
  
PACKAGE_ADD_EXECUTABLE(
  SaPtest
  SOURCES SaPtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )