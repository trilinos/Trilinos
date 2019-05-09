TRIBITS_COPY_FILES_TO_BINARY_DIR(YAML_TEST_FILES
  SOURCE_FILES "plasma_oscillation_rtc.xml.yaml"
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Parser_UnitTests  
  SOURCES
    unit_tests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
