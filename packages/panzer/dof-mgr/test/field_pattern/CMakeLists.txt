INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

IF (Intrepid2_ENABLE_KokkosDynRankView)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tGeometricAggFieldPattern.cpp
  SOURCES tGeometricAggFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFieldPattern
  SOURCES tFieldPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFieldAggPattern
  SOURCES tFieldAggPattern.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

ENDIF()
