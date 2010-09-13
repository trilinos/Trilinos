
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tFieldPattern
  SOURCES tFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tGeometricAggFieldPattern
  SOURCES tGeometricAggFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tIntrepidFieldPattern
  SOURCES tIntrepidFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tFieldAggPattern
  SOURCES tFieldAggPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )
