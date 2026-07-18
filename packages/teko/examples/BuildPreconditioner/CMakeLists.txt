ADD_SUBDIRECTORY(step1)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  BuildPreconditioner_Driver_Example
  SOURCES example-driver.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  BuildPreconditioner_Driver_Example-Belos
  SOURCES example-driver-belos.cpp
  COMM serial mpi
)
