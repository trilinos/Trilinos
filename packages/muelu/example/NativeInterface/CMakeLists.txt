INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  PlainAggregation
  SOURCES PlainAggregation.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
)
	