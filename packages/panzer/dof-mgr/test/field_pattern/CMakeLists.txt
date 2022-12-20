TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tGeometricAggFieldPattern2
  SOURCES tGeometricAggFieldPattern2.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFieldPattern2
  SOURCES tFieldPattern2.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFieldAggPattern2
  SOURCES tFieldAggPattern2.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFieldAggPattern_DG
  SOURCES tFieldAggPattern_DG.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEntityFieldPattern
  SOURCES tEntityFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )
