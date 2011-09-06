INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageAddAdvancedTest)

ADD_DEFINITIONS(-DUSING_CMAKE)

find_program (${PACKAGE_NAME}_Netcdf_EXECUTABLE_NCDUMP
	     NAME ncdump
	     PATHS 
	        $ENV{ACCESS}/bin
	        $ENV{NETCDF_DIR}/bin
	        $ENV{NETCDF_DIR}/ncdump
		ENV PATH
	     )
	     
# ===============================================


IF (${PACKAGE_NAME}_Netcdf_EXECUTABLE_NCDUMP)
  SET( forbind_OLDTESTS
    testwt
    testwtd
    testrd
    testrd1
    testrdd
    testwt1
    testwt2
    testwtm
    testcp
    testcpnl
    testcpln
    testwt_nsid
    testrd_nsid
  )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
    ${CMAKE_CURRENT_BINARY_DIR}/testall
    @ONLY
  )

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved depedencies
  #       when built with shareable libraries on some platforms....

  FOREACH ( testName ${forbind_OLDTESTS} )
    PACKAGE_ADD_EXECUTABLE( f_${testName} NOEXEPREFIX LINKER_LANGUAGE Fortran SOURCES ${testName}.f DEPLIBS exodus_for)
  ENDFOREACH ( testName )

  PACKAGE_ADD_ADVANCED_TEST(
   exodus_unit_tests
   TEST_0 CMND /bin/csh ARGS ${CMAKE_CURRENT_BINARY_DIR}/testall
	  OUTPUT_FILE exodus_unit_tests.out
	  PASS_ANY
   TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_unit_tests.gold 
		         ${CMAKE_CURRENT_BINARY_DIR}/exodus_unit_tests.out
   COMM serial		      
   FINAL_PASS_REGULAR_EXPRESSION 
  )
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus unit tests will not be run due to failure finding ncdump executable.")
ENDIF()
