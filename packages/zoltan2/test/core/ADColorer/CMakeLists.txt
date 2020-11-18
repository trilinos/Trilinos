
TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ADColorer
      NAME ADColorer_simple
      SOURCES ADColorer.cpp
      ARGS "--file=simple"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

IF (${PROJECT_NAME}_ENABLE_Galeri)
TRIBITS_ADD_TEST(
    ADColorer
    NAME ADColorer_galeri
    ARGS "--xdim=10 --ydim=8 --zdim=6"
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)
ENDIF()
  
# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_adcolor_tests
  SOURCE_FILES
    simple.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

