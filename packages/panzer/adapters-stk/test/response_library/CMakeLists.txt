

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/disc-fe/test/equation_set)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/disc-fe/test/closure_model)
INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/disc-fe/test/response_library)
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  STK_ResponseLibraryTest2
  SOURCES response_library2.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  STK_VolumeSideResponse
  SOURCES volumetric_side_response.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )
