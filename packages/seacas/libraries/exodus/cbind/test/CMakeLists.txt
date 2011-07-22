INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageAddAdvancedTest)

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

PACKAGE_ADD_EXECUTABLE(
 ExoIICTests
 SOURCES ${ExoIICTests}
)

find_program (${PACKAGE_NAME}_Netcdf_EXECUTABLE_NCDUMP
	     NAME ncdump
	     PATHS 
	        $ENV{ACCESS}/bin
	        $ENV{NETCDF_DIR}/bin
	        $ENV{NETCDF_DIR}/ncdump
		ENV PATH
	     )
	     
PACKAGE_ADD_TEST(
	ExoIICTests
	NAME CreateEdgeFaceWithConcats
	COMM mpi serial
	NUM_MPI_PROCS 1
	ARGS "CreateEdgeFace -pcab -pcset -pvpax"
)

PACKAGE_ADD_TEST(
	ExoIICTests
	NAME ReadEdgeFaceWithConcats
	COMM mpi serial
	NUM_MPI_PROCS 1
	ARGS ReadEdgeFace
)

# Should be a better way to do this, but...
if (TPL_ENABLE_MPI)
 set_property(TEST SEACAS_ReadEdgeFaceWithConcats_MPI_1 APPEND PROPERTY DEPENDS SEACAS_CreateEdgeFaceWithConcats_MPI_1)
ELSE()
 set_property(TEST SEACAS_ReadEdgeFaceWithConcats APPEND PROPERTY DEPENDS SEACAS_CreateEdgeFaceWithConcats)
ENDIF()

# ===============================================


IF (${PACKAGE_NAME}_Netcdf_EXECUTABLE_NCDUMP)
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
  testcp_ln
  testcp_nl
  testcpd
  testwt_clb
  testwt_nc
  testrd_nc
  testwt-zeroe
  testwt-zeron
  testrdv
  testwt-long-name
  testrd-long-name
  testwt-one-attrib
  testwt-partial
  testrd-partial
  testrd-nsided
  testwt-nsided
  testrd-nfaced
  testwt-nfaced
  create_mesh
  rd_wt_mesh
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

FOREACH ( testName ${cbind_OLDTESTS} )
  PACKAGE_ADD_EXECUTABLE( ${testName} NOEXEPREFIX SOURCES ${testName}.c LINKER_LANGUAGE CXX DEPLIBS exodus)
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
