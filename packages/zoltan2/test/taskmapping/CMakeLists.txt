
#IF (${PROJECT_NAME}_ENABLE_SEACASExodus)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
      taskMappingFromExodus
      SOURCES taskMappingFromExodus.cpp
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

#Use the macro below if need to copy files for testing
#TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_taskmapping_tests
#  SOURCE_FILES
#    filename1.exoII
#    filename2.exoII
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#)

#ENDIF()
