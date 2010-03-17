
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../stk_util)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/src)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
UnitTestBucket.cpp
UnitTestBoundaryAnalysis.cpp
UnitTestBoxMeshFixture.cpp
UnitTestBulkData_ChangeParts.cpp
UnitTestBulkData.cpp
UnitTestBulkData_Destroy.cpp
UnitTestBulkData_Generate.cpp
UnitTestBulkData_new.cpp
UnitTestBulkModification.cpp
UnitTestDataTraits.cpp
UnitTestField.cpp
UnitTestFixture.cpp
UnitTestGetBuckets.cpp
UnitTestGridMeshFixture.cpp
UnitTestMain.cpp
UnitTestMesh.cpp
UnitTestMetaData.cpp
UnitTestRelation.cpp
UnitTestRingMeshFixture.cpp
UnitTestSelector.cpp
UnitTestTemplate.cpp
UnitTestTransactionLog.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    stk_mesh_unit_tests
    SOURCES ${SOURCES}
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

