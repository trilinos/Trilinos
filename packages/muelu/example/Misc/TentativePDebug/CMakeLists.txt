INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../adapters/belos)
#######################################################################################
#######################################################################################

TRIBITS_ADD_EXECUTABLE(
  tentativePtest
  SOURCES tentativePtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )

TRIBITS_ADD_EXECUTABLE(
  SaPtest
  SOURCES SaPtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )

TRIBITS_ADD_EXECUTABLE(
  PgPtest
  SOURCES PgPtest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )
