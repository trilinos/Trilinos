INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tFullRelations
  SOURCES tFullRelations.cpp Panzer_Relations.cpp Panzer_FaceFieldPattern.cpp Panzer_EdgeFieldPattern.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

