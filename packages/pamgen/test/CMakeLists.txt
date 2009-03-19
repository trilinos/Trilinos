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
# Set up the test arguments
#

SET( NUM_TESTS 2 )

SET( CMND_ARGS_ARRAY
  "-f mesh.txt -r 1 -n 3 -d 3"
  "-f mesh.txt -r 2 -n 3 -d 3"
  )


#
# Set up the tests
#

SET( RES_FILES )

FOREACH( TEST_IDX RANGE 1 ${NUM_TESTS} )

  SET( FILE_IDX_EXT 0${TEST_IDX} )  # ToDo: Generalize for above 09!

  SET( RES_FILE run_test_script_${FILE_IDX_EXT}.res )

  LIST( APPEND RES_FILES ${RES_FILE} )

  MATH( EXPR ARRAY_IDX ${TEST_IDX}-1 )

  LIST( GET CMND_ARGS_ARRAY ${ARRAY_IDX} CMND_ARGS )

  SET( TMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/temp_file_${FILE_IDX_EXT} )

  PACKAGE_ADD_ADVANCED_TEST(
    run_test_script_${FILE_IDX_EXT}
    TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/pamgen_test_driver.exe
      ARGS ${CMND_ARGS}
      OUTPUT_FILE ${TMP_FILE}
    TEST_1 CMND diff
      ARGS ${TMP_FILE} ${RES_FILE}
    )

ENDFOREACH()


PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES mesh.txt ${RES_FILES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
  )
