INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)


PACKAGE_ADD_EXECUTABLE(
  pamgen_test_driver
  NOEXEPREFIX
  SOURCES
    ../example/pamgen_test_driver.c
    ../example/getopts.C
  )


#
# Set up a simple function that will add the tests
#

SET( RES_FILES )


FUNCTION(PAMGEN_ADD_TEST  FILE_IDX_EXT)

  SET( RES_FILE run_test_script_${FILE_IDX_EXT}.res )

  SET( RES_FILES ${RES_FILES} ${RES_FILE} PARENT_SCOPE )

  SET( TMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/temp_file_${FILE_IDX_EXT} )

  PACKAGE_ADD_ADVANCED_TEST(
    run_test_script_${FILE_IDX_EXT}
    TEST_0 EXEC pamgen_test_driver NOEXEPREFIX
      ARGS ${ARGN}
      OUTPUT_FILE ${TMP_FILE}
      NO_ECHO_OUTPUT
    TEST_1 CMND diff
      ARGS ${TMP_FILE} ${RES_FILE}
    COMM serial
    )

ENDFUNCTION()


#
# Define the tests
#

PAMGEN_ADD_TEST( 01 -f mesh.txt -r 1 -n 3 -d 3 )
PAMGEN_ADD_TEST( 02 -f mesh.txt -r 2 -n 3 -d 3 )


#
# Copy the files over for each test
#

PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES mesh.txt ${RES_FILES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
  )
