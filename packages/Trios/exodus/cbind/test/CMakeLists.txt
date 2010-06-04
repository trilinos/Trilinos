INCLUDE(PackageAddExecutableAndTest)

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
 COMM serial
)

PACKAGE_ADD_TEST(
	ExoIICTests
	NAME CreateEdgeFaceWithConcats
	NUM_MPI_PROCS 1
	COMM serial
	ARGS "CreateEdgeFace -pcab -pcset -pvpax"
)

PACKAGE_ADD_TEST(
	ExoIICTests
	NAME ReadEdgeFaceWithConcats
	NUM_MPI_PROCS 1
	COMM serial
	ARGS ReadEdgeFace
)

# ===============================================

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

FOREACH ( testName ${cbind_OLDTESTS} )
  ADD_EXECUTABLE( ${testName} ${testName}.c )
  TARGET_LINK_LIBRARIES( ${testName} exodus )
ENDFOREACH ( testName )

ADD_CUSTOM_TARGET( check
  COMMAND /bin/csh ${CMAKE_CURRENT_BINARY_DIR}/testall
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${cbind_OLDTESTS}
)

