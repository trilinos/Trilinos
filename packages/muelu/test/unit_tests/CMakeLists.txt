
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
main.cpp
Level.cpp
Hierarchy.cpp
RAPFactory.cpp
SaPFactory.cpp
TransPFactory.cpp
SmootherFactory.cpp
IfpackSmoother.cpp
Needs.cpp
GenericPRFactory.cpp
Amesos.cpp
TentativePFactory.cpp
UCAggregationFactory.cpp
)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial  mpi
    DEPLIBS muelu muelu-cthulhu
  )
