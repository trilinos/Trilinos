INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../adapters/belos)
#######################################################################################
#######################################################################################

PACKAGE_ADD_EXECUTABLE(
  tentativePtest
  SOURCES tentativePtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )
  
PACKAGE_ADD_EXECUTABLE(
  SaPtest
  SOURCES SaPtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )
  
PACKAGE_ADD_EXECUTABLE(
  PgPtest
  SOURCES PgPtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )