INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddAdvancedTest)
INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(ParseVariableArguments)


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

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "HOST;XHOST;HOSTTYPE;XHOSTTYPE;ARGS"
     #options
     "LONG_LONG"
     ${ARGN}
     )

  SET(HOST_ARGS "")
  IF(PARSE_HOST)
    SET(HOST_ARGS ${HOST_ARGS} HOST ${PARSE_HOST})
  ENDIF()
  IF(PARSE_XHOST)
    SET(HOST_ARGS ${HOST_ARGS} XHOST ${PARSE_XHOST})
  ENDIF()
  IF(PARSE_HOSTTYPE)
    SET(HOST_ARGS ${HOST_ARGS} HOSTTYPE ${PARSE_HOSTTYPE})
  ENDIF()
  IF(PARSE_XHOSTTYPE)
    SET(HOST_ARGS ${HOST_ARGS} XHOSTTYPE ${PARSE_XHOSTTYPE})
  ENDIF()
  
  #if LONG_LONG is set then we need to use the pamgen_test_driver_l driver to run tests
  #this generally means just adding "_l" to a few files
  SET(INT_TYPE "")
  IF(PARSE_LONG_LONG)
    SET(INT_TYPE "_l")
  ENDIF()

  SET( RES_FILE run_test_script${INT_TYPE}_${FILE_IDX_EXT}.res )

  SET( RES_FILES ${RES_FILES} ${RES_FILE} PARENT_SCOPE )

  SET( TMP_FILE ${CMAKE_CURRENT_BINARY_DIR}/temp_file${INT_TYPE}_${FILE_IDX_EXT} )

  PACKAGE_ADD_ADVANCED_TEST(
    run_test_script${INT_TYPE}_${FILE_IDX_EXT}
    TEST_0 EXEC pamgen_test_driver${INT_TYPE} NOEXEPREFIX
      ARGS ${PARSE_ARGS}
      OUTPUT_FILE ${TMP_FILE}
      NO_ECHO_OUTPUT
    TEST_1 CMND diff
      ARGS ${TMP_FILE} ${RES_FILE}
    COMM serial
    ${HOST_ARGS}
    )

ENDFUNCTION()

#
# Define the tests
#

ASSERT_DEFINED(CMAKE_BUILD_TYPE)
IF (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  PAMGEN_ADD_TEST(01 LONG_LONG ARGS -f lbrick.txt -r 9999 -n 10000 -d 3 )
ENDIF()


PAMGEN_ADD_TEST( 01 ARGS -f mesh.txt -r 1 -n 3 -d 3 )
PAMGEN_ADD_TEST( 02 ARGS -f mesh.txt -r 2 -n 3 -d 3 )
PAMGEN_ADD_TEST( 03 ARGS -f mesh.txt -r 0 -n 3 -d 3 )
PAMGEN_ADD_TEST( 04 ARGS -f brick.txt -r 1 -n 3 -d 3 )
PAMGEN_ADD_TEST( 05 ARGS -f brick.txt -r 2 -n 3 -d 3 )
PAMGEN_ADD_TEST( 06 ARGS -f brick.txt -r 0 -n 3 -d 3 )
PAMGEN_ADD_TEST( 07 ARGS -f tri.txt -r 0 -n 6 -d 3 )
PAMGEN_ADD_TEST( 08 ARGS -f tri.txt -r 5 -n 6 -d 3 )
PAMGEN_ADD_TEST( 09 ARGS -f tri_proc.txt -r 5 -n 6 -d 3 )
PAMGEN_ADD_TEST( 10 ARGS -f 2dtri.txt -r 0 -n 1 -d 2 )
PAMGEN_ADD_TEST( 11 ARGS -f 2dtri.txt -r 0 -n 2 -d 2 )
PAMGEN_ADD_TEST( 13 ARGS -f 2dtri.txt -r 3 -n 4 -d 2 )
PAMGEN_ADD_TEST( 14 ARGS -f 2dtribis.txt -r 0 -n 2 -d 2 )
PAMGEN_ADD_TEST( 15 ARGS -f 2dtribis.txt -r 3 -n 4 -d 2 )
PAMGEN_ADD_TEST( 16 ARGS -f hbrick.txt -r 0 -n 100 -d 3 XHOSTTYPE AIX )
IF (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  PAMGEN_ADD_TEST( 17 ARGS -f 3dbj.txt -r 0 -n 1 -d 3 )
ENDIF()
PAMGEN_ADD_TEST( 18 ARGS -f hbrick.txt -r 11 -n 100 -d 3 XHOSTTYPE AIX )
PAMGEN_ADD_TEST( 19 ARGS -f trisect_multi_theta.txt )
PAMGEN_ADD_TEST( 20 ARGS -f trisect_multi_theta_2d.txt -d 2)
PAMGEN_ADD_TEST( 21 ARGS -f trisect_multi_theta_2d.txt --all -n 6 -d 2)
PAMGEN_ADD_TEST( 22 ARGS -f ui.txt --all -n 4 )
PAMGEN_ADD_TEST( 23 ARGS -f ui2d.txt --all -n 4 -d 2)
PAMGEN_ADD_TEST( 24 ARGS -f geom_warp.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 25 ARGS -f geom_warp2d.txt -r 0 -n 1 -d 2)
IF (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  PAMGEN_ADD_TEST( 26 ARGS -f report1.txt -r 0 -n 1 -d 3)
ENDIF()
PAMGEN_ADD_TEST( 27 ARGS -f report2.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 28 ARGS -f report3.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 29 ARGS -f report4.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 30 ARGS -f report5.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 31 ARGS -f report6.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 32 ARGS -f report7.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 33 ARGS -f report8.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 34 ARGS -f report9.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 35 ARGS -f report10.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 36 ARGS -f report11.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 37 ARGS -f report12.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 38 ARGS -f report13.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 39 ARGS -f report14.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 40 ARGS -f report15.txt -r 0 -n 1 -d 2)
PAMGEN_ADD_TEST( 41 ARGS -f report16.txt -r 0 -n 8 -d 3)
PAMGEN_ADD_TEST( 42 ARGS -f report17.txt -r 0 -n 8 -d 3)
PAMGEN_ADD_TEST( 43 ARGS -f uded.txt -r 0 -n 1 -d 3)
PAMGEN_ADD_TEST( 44 ARGS -f ur.txt -r 0 -n 1 -d 3)

# 2009/05/28: rabartl: Above, the test 17 timing out at 600s on gabriel and
# the test 26 just takes too long (see bug 4428).


#
# Copy the files over for each test
#

PACKAGE_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  SOURCE_FILES geom_warp2d.txt uded.txt ur.txt  report17.txt report16.txt report15.txt report14.txt report13.txt report12.txt report11.txt report9.txt report10.txt report7.txt report8.txt report6.txt report4.txt report5.txt report3.txt report2.txt report1.txt geom_warp.txt ui2d.txt ui.txt mesh.txt lbrick.txt brick.txt tri.txt tri_proc.txt 2dtri.txt 2dtribis.txt hbrick.txt 3dbj.txt trisect_multi_theta.txt trisect_multi_theta_2d.txt ${RES_FILES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
  )
