TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/unit_tests)

TRIBITS_ADD_EXECUTABLE(
  RegionsDriver
  SOURCES Driver.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  StructuredRegionsDriver
  SOURCES Driver_Structured_Regions.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(tawiesn_cp SOURCE_FILES structured_1dof.xml)
