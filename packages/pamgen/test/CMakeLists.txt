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

PACKAGE_ADD_EXECUTABLE(
  pamgen_test_driver_l
  NOEXEPREFIX
  SOURCES
    ../example/pamgen_test_driver_l.c
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

FUNCTION(PAMGEN_ADD_TEST_L  FILE_IDX_EXT)

  SET( RES_FILE run_test_script_l_${FILE_IDX_EXT}.res )

  SET( RES_FILES ${RES_FILES} ${RES_FILE} PARENT_SCOPE )

  SET( TMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/temp_file_l_${FILE_IDX_EXT} )

  PACKAGE_ADD_ADVANCED_TEST(
    run_test_script_l_${FILE_IDX_EXT}
    TEST_0 EXEC pamgen_test_driver_l NOEXEPREFIX
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

ASSERT_DEFINED(CMAKE_BUILD_TYPE)

PAMGEN_ADD_TEST_L( 01 -f lbrick.txt -r 9999 -n 10000 -d 3 )

PAMGEN_ADD_TEST( 01 -f mesh.txt -r 1 -n 3 -d 3 )
PAMGEN_ADD_TEST( 02 -f mesh.txt -r 2 -n 3 -d 3 )
PAMGEN_ADD_TEST( 03 -f mesh.txt -r 0 -n 3 -d 3 )
PAMGEN_ADD_TEST( 04 -f brick.txt -r 1 -n 3 -d 3 )
PAMGEN_ADD_TEST( 05 -f brick.txt -r 2 -n 3 -d 3 )
PAMGEN_ADD_TEST( 06 -f brick.txt -r 0 -n 3 -d 3 )
PAMGEN_ADD_TEST( 07 -f tri.txt -r 0 -n 6 -d 3 )
PAMGEN_ADD_TEST( 08 -f tri.txt -r 5 -n 6 -d 3 )
PAMGEN_ADD_TEST( 09 -f tri_proc.txt -r 5 -n 6 -d 3 )
PAMGEN_ADD_TEST( 10 -f 2dtri.txt -r 0 -n 1 -d 2 )
PAMGEN_ADD_TEST( 11 -f 2dtri.txt -r 0 -n 2 -d 2 )
PAMGEN_ADD_TEST( 13 -f 2dtri.txt -r 3 -n 4 -d 2 )
PAMGEN_ADD_TEST( 14 -f 2dtribis.txt -r 0 -n 2 -d 2 )
PAMGEN_ADD_TEST( 15 -f 2dtribis.txt -r 3 -n 4 -d 2 )
PAMGEN_ADD_TEST( 16 -f hbrick.txt -r 0 -n 100 -d 3 )
IF (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  PAMGEN_ADD_TEST( 17 -f 3dbj.txt -r 0 -n 1 -d 3 )
ENDIF()
PAMGEN_ADD_TEST( 18 -f hbrick.txt -r 11 -n 100 -d 3 )
PAMGEN_ADD_TEST( 19 -f trisect_multi_theta.txt )
PAMGEN_ADD_TEST( 20 -f trisect_multi_theta_2d.txt -d 2)
PAMGEN_ADD_TEST( 21 -f trisect_multi_theta_2d.txt --all -n 6 -d 2)
PAMGEN_ADD_TEST( 22 -f ui.txt --all -n 4 )
PAMGEN_ADD_TEST( 23 -f ui2d.txt --all -n 4 -d 2)
PAMGEN_ADD_TEST( 24 -f geom_warp.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 25 -f geom_warp2d.txt -r 0 -n 1 -d 2)
IF (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  PAMGEN_ADD_TEST( 26 -f report1.txt -r 0 -n 1 -d 3)
ENDIF()
PAMGEN_ADD_TEST( 27 -f report2.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 28 -f report3.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 29 -f report4.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 30 -f report5.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 31 -f report6.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 32 -f report7.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 33 -f report8.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 34 -f report9.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 35 -f report10.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 36 -f report11.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 37 -f report12.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 38 -f report13.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 39 -f report14.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 40 -f report15.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 41 -f report16.txt -r 0 -n 8 -d 3)
PAMGEN_ADD_TEST( 42 -f report17.txt -r 0 -n 8 -d 3)
PAMGEN_ADD_TEST( 43 -f uded.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 44 -f ur.txt -r 0 -n 1 -d 3)

# 2009/05/28: rabartl: Above, the test 17 timing out at 600s on gabriel and
# the test 26 just takes too long (see bug 4428).


#
# Copy the files over for each test
#

PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES geom_warp2d.txt uded.txt ur.txt  report17.txt report16.txt report15.txt report14.txt report13.txt report12.txt report11.txt report9.txt report10.txt report7.txt report8.txt report6.txt report4.txt report5.txt report3.txt report2.txt report1.txt geom_warp.txt ui2d.txt ui.txt mesh.txt lbrick.txt brick.txt tri.txt tri_proc.txt 2dtri.txt 2dtribis.txt hbrick.txt 3dbj.txt trisect_multi_theta.txt trisect_multi_theta_2d.txt ${RES_FILES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
  )
