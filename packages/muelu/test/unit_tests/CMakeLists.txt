
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

#PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

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
  Needs.cpp
  GenericPRFactory.cpp
  TentativePFactory.cpp
  Aggregates.cpp
  CoalesceDropFactory.cpp
  Memory.cpp
  #UCAggregationFactory.cpp
  ExtendedHashtable.cpp
)

### Tests that require other Trilinos packages

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Ifpack
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  )

IF (${PACKAGE_NAME}_ENABLE_Ifpack2)
  APPEND_SET(SOURCES
  Ifpack2Smoother.cpp
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(SOURCES
  IfpackSmoother.cpp
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(SOURCES
  Amesos.cpp
  )
ENDIF()

PACKAGE_ADD_EXECUTABLE(
  UnitTests
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
  )

PACKAGE_ADD_TEST(
  UnitTests
  NAME "UnitTestsEpetra"
  ARGS "--linAlgebra=0"
  PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

PACKAGE_ADD_TEST(
  UnitTests
  NAME "UnitTestsTpetra"
  ARGS "--linAlgebra=1"
  PASS_REGULAR_EXPRESSION "End Result: TEST PASSED"
  NUM_MPI_PROCS 4
  COMM serial mpi
  )
