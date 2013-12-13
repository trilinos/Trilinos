INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  test_version
  SOURCES version.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(core_tests
#  SOURCE_FILES parameter_list_acceptance_test.xml parameter_list_acceptance_test2.xml
#  EXEDEPS parameter_list_acceptance_test
#  )
