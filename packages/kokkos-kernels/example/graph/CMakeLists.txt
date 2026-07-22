include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES "PartitioningExample.cpp")

tribits_add_executable(PartitioningExample SOURCES ${SOURCES} COMM serial)
