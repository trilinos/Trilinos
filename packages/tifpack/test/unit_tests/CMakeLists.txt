
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/NodeAPI)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/LinAlg)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../kokkos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../tpetra/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../tpetra/src)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
Tifpack_UnitTestMain.cpp
Tifpack_UnitTestTemplate.cpp
Tifpack_UnitTestCreateOverlapGraph.cpp
Tifpack_UnitTestOverlapGraph.cpp
Tifpack_UnitTestParameters.cpp
Tifpack_UnitTestIlukGraph.cpp
Tifpack_UnitTestILUT.cpp
Tifpack_UnitTestPointRelaxation.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tifpack_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    DEPLIBS tpetra kokkos teuchos
  )

