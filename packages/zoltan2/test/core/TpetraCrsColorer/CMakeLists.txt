
TRIBITS_ADD_EXECUTABLE_AND_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple
      SOURCES TpetraCrsColorer.cpp
      ARGS "--file=simple"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

TRIBITS_ADD_TEST(
      TpetraCrsColorer
      NAME TpetraCrsColorer_simple_nodistrib
      ARGS "--file=simple --no-distribute"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

IF (${PROJECT_NAME}_ENABLE_Galeri)
TRIBITS_ADD_TEST(
    TpetraCrsColorer
    NAME TpetraCrsColorer_galeri
    ARGS "--xdim=20 --ydim=16 --zdim=10"
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

