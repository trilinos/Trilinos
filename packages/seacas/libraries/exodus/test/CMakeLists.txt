
ADD_DEFINITIONS(-DUSING_CMAKE)

SET( cbind_TEST_SRCS
  CreateEdgeFace.c
  ReadEdgeFace.c
)

# Always report error/warning messages during testing:
# The ex_opts call precipitates the EXTRA_INCLUDE below.
SET( CMAKE_TESTDRIVER_BEFORE_TESTMAIN
  "  ex_opts( EX_VERBOSE );\n"
)
CREATE_TEST_SOURCELIST( ExoIICTests
  ExoIICTests.cxx
  ${cbind_TEST_SRCS}
  EXTRA_INCLUDE exodusII_test.h
  )

TRIBITS_ADD_EXECUTABLE(
  ExoIICTests
  SOURCES ${ExoIICTests}
)

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
	     
TRIBITS_ADD_TEST(
	ExoIICTests
	NAME CreateEdgeFaceWithConcats
	COMM mpi serial
	NUM_MPI_PROCS 1
	ARGS "CreateEdgeFace -pcab -pcset -pvpax"
)

TRIBITS_ADD_TEST(
	ExoIICTests
	NAME ReadEdgeFaceWithConcats
	COMM mpi serial
	NUM_MPI_PROCS 1
	ARGS ReadEdgeFace
)

# Should be a better way to do this, but...
if (TPL_ENABLE_MPI)
  set_property(TEST ${PACKAGE_NAME}_ReadEdgeFaceWithConcats_MPI_1 APPEND PROPERTY DEPENDS ${PACKAGE_NAME}_CreateEdgeFaceWithConcats_MPI_1)
ELSE()
  set_property(TEST ${PACKAGE_NAME}_ReadEdgeFaceWithConcats APPEND PROPERTY DEPENDS ${PACKAGE_NAME}_CreateEdgeFaceWithConcats)
ENDIF()

# ===============================================

IF ( NETCDF_NCDUMP_BINARY )
  SET(NETCDF_NCDUMP ${NETCDF_NCDUMP_BINARY})
  SET( cbind_OLDTESTS
    testwt
    testwtd
    testrd
    testrd1
    testrdd
    testwt1
    testwt2
    testwtm
    testwt_ss
    testwt_nossnsdf
    testrd_ss
    testrdwt
    testcp
    testcp_nl
    testcpd
    testwt_clb
    testwt_nc
    testrd_nc
    testwt-zeroe
    testwt-zeron
    testwt-long-name
    testrd-long-name
    testwt-one-attrib
    testwt-partial
    testrd-nsided
    testwt-nsided
    testrd-nfaced
    testwt-nfaced
    test_nemesis
    create_mesh
    rd_wt_mesh
    test-empty
    testwt-compress
  )

  IF (SEACASExodus_ENABLE_THREADSAFE)
     SET(TS "_ts")
     SET(cbind_OLDTESTS
       ${cbind_OLDTESTS}
       test_ts_nvar
       test_ts_nvar_rd
       test_ts_partial_nvar
       test_ts_partial_nvar_rd
       test_ts_files
       test_ts_errval
     )
  ENDIF()

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
    ${CMAKE_CURRENT_BINARY_DIR}/testall
    @ONLY
  )

  if ( TPL_Netcdf_Enables_Netcdf4 )
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testnc4.in
    ${CMAKE_CURRENT_BINARY_DIR}/testnc4
    @ONLY
  )
  endif()

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved depedencies
  #       when built with shareable libraries on some platforms....

  FOREACH ( testName ${cbind_OLDTESTS} )
    TRIBITS_ADD_EXECUTABLE( ${testName} NOEXEPREFIX NOEXESUFFIX SOURCES ${testName}.c LINKER_LANGUAGE CXX)
  ENDFOREACH ( testName )

  TRIBITS_ADD_ADVANCED_TEST(
   exodus_unit_tests
   TEST_0 CMND /bin/bash ARGS ${CMAKE_CURRENT_BINARY_DIR}/testall 
	  OUTPUT_FILE exodus_unit_tests.out
	  PASS_ANY
   TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_unit_tests.gold${TS}
		         ${CMAKE_CURRENT_BINARY_DIR}/exodus_unit_tests.out
   COMM mpi serial
   OVERALL_NUM_MPI_PROCS 1
   FINAL_PASS_REGULAR_EXPRESSION 
  )

  if ( TPL_Netcdf_Enables_Netcdf4 )
    TRIBITS_ADD_ADVANCED_TEST(
     exodus_unit_tests_nc4_env
     TEST_0 CMND /bin/bash ARGS ${CMAKE_CURRENT_BINARY_DIR}/testall netcdf4 
	    OUTPUT_FILE exodus_unit_tests.out
	    PASS_ANY
     TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_unit_tests.gold${TS} 
                           ${CMAKE_CURRENT_BINARY_DIR}/exodus_unit_tests.out
     COMM mpi serial
     FINAL_PASS_REGULAR_EXPRESSION 
    )

    TRIBITS_ADD_ADVANCED_TEST(
     exodus_nc4_unit_tests
     TEST_0 CMND /bin/bash ARGS ${CMAKE_CURRENT_BINARY_DIR}/testnc4
	    OUTPUT_FILE exodus_nc4_unit_tests.out
	    PASS_ANY
     TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_nc4_unit_tests.gold 
		           ${CMAKE_CURRENT_BINARY_DIR}/exodus_nc4_unit_tests.out
     COMM mpi serial
     OVERALL_NUM_MPI_PROCS 1
     FINAL_PASS_REGULAR_EXPRESSION 
    )
  endif()

#  if ( TPL_Netcdf_Enables_Netcdf5 )  
    TRIBITS_ADD_ADVANCED_TEST(
     exodus_unit_tests_nc5_env
     TEST_0 CMND /bin/bash ARGS ${CMAKE_CURRENT_BINARY_DIR}/testall netcdf5
	    OUTPUT_FILE exodus_unit_tests.out
	    PASS_ANY
     TEST_1 CMND diff ARGS ${CMAKE_CURRENT_SOURCE_DIR}/exodus_unit_tests.gold${TS} 
                           ${CMAKE_CURRENT_BINARY_DIR}/exodus_unit_tests.out
     COMM mpi serial
     FINAL_PASS_REGULAR_EXPRESSION 
    )
#  endif()
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus unit tests will not be run due to failure finding ncdump executable.")
ENDIF()
