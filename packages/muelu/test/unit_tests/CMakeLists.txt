
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
  MueLu_UnitTests.cpp
  MueLu_TestHelpers.cpp
  Level.cpp
  Hierarchy.cpp
  RAPFactory.cpp
  SaPFactory.cpp
  TransPFactory.cpp
  SmootherFactory.cpp
  IfpackSmoother.cpp
  Ifpack2Smoother.cpp
  Needs.cpp
  GenericPRFactory.cpp
  Amesos.cpp
  TentativePFactory.cpp
  Aggregates.cpp
  UCAggregationFactory.cpp
  CoalesceDropFactory.cpp
)

PACKAGE_ADD_EXECUTABLE(
  UnitTests
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS muelu muelu-cthulhu
  )

PACKAGE_ADD_TEST(
  UnitTests
  NAME "UnitTestsEpetra"
  ARGS "--linAlgebra=0"
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

PACKAGE_ADD_TEST(
  UnitTests
  NAME "UnitTestsTpetra"
  ARGS "--linAlgebra=1"
  NUM_MPI_PROCS 4
  COMM serial mpi
  )