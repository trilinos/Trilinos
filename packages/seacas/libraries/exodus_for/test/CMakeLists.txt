ADD_DEFINITIONS(-DUSING_CMAKE)

if ( NETCDF_NCDUMP_BINARY )
   # Do nothing. Set earlier in config process...
else()
   find_program (NETCDF_NCDUMP_BINARY
	     NAME ncdump
	     PATHS
	       ${NetCDF_BINARY_DIR}
               ${Netcdf_LIBRARY_DIRS}/../bin
               $ENV{ACCESS}/bin
               $ENV{NETCDF_DIR}/bin
               $ENV{NETCDF_DIR}/ncdump
               ENV PATH
	        ${NetCDF_DIR}/bin
	     )
endif()

# ===============================================

IF ( NETCDF_NCDUMP_BINARY )
  SET(NETCDF_NCDUMP ${NETCDF_NCDUMP_BINARY})
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
    testwt_nsid
    testrd_nsid
    test_nem
  )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
    ${CMAKE_CURRENT_BINARY_DIR}/testall
    @ONLY
  )

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved dependencies
  #       when built with shareable libraries on some platforms....

  FOREACH ( testName ${forbind_OLDTESTS} )
    TRIBITS_ADD_EXECUTABLE( f_${testName} NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES ${testName}.f)
  ENDFOREACH ( testName )

  TRIBITS_ADD_ADVANCED_TEST(
   exodus_unit_tests
   TEST_0 CMND /bin/bash ARGS ${CMAKE_CURRENT_BINARY_DIR}/testall
	  OUTPUT_FILE exodus_unit_tests.out
	  PASS_ANY
   TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_unit_tests.gold
		         ${CMAKE_CURRENT_BINARY_DIR}/exodus_unit_tests.out
   COMM mpi serial
   OVERALL_NUM_MPI_PROCS 1
   FINAL_PASS_REGULAR_EXPRESSION
  )
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus unit tests will not be run due to failure finding ncdump executable.")
ENDIF()
