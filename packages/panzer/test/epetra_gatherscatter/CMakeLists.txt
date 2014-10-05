
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetraGather
  SOURCES tEpetraGather.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetraScatter
  SOURCES tEpetraScatter.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetraScatterDirichlet
  SOURCES tEpetraScatterDirichlet.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )
