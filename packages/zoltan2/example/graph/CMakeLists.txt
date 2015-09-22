
IF (${PACKAGE_NAME}_ENABLE_Scotch AND ${PACKAGE_NAME}_ENABLE_ParMETIS)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  graph_example
  NAME parmetis_example
  SOURCES graph.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS "--method=parmetis"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  graph_example
  NAME scotch_example
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS "--method=scotch"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()
