
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
        ADDED_TESTS_NAMES_OUT CreateEdgeFaceWithConcatsAddedName
)

TRIBITS_ADD_TEST(
	ExoIICTests
	NAME ReadEdgeFaceWithConcats
	COMM mpi serial
	NUM_MPI_PROCS 1
	ARGS ReadEdgeFace
        ADDED_TESTS_NAMES_OUT ReadEdgeFaceWithConcatsAddedName
)

if (CreateEdgeFaceWithConcatsAddedName AND ReadEdgeFaceWithConcatsAddedName)
  set_property(TEST ${ReadEdgeFaceWithConcatsAddedName}
   APPEND PROPERTY DEPENDS ${CreateEdgeFaceWithConcatsAddedName})
endif()

# ===============================================

IF ( NETCDF_NCDUMP_BINARY )
  SET(NETCDF_NCDUMP ${NETCDF_NCDUMP_BINARY})
  SET( cbind_OLDTESTS
    test-field-utils
    testwt
    testwtd
    testwtd-to-f
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
    testcp_tran
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
    testwt-results
    testwt-oned
    testwt-assembly
    testrd-assembly
    test-add-assembly
    testwt-blob
    testrd-blob
    testwt-field-metadata
    testrd-field-metadata
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
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall
    @ONLY
  )
  # now copy the temporary into the final destination, setting the permissions
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  if ( TPL_Netcdf_Enables_Netcdf4 )
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testnc4.in
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testnc4
    @ONLY
  )
  # now copy the temporary into the final destination, setting the permissions
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testnc4
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved dependencies
  #       when built with shareable libraries on some platforms....

  FOREACH ( testName ${cbind_OLDTESTS} )
    TRIBITS_ADD_EXECUTABLE( ${testName} NOEXEPREFIX NOEXESUFFIX SOURCES ${testName}.c LINKER_LANGUAGE CXX)
    set_property(TARGET ${testName} PROPERTY C_STANDARD 99)
  ENDFOREACH ( testName )

  TRIBITS_ADD_TEST(
     test-field-utils
     NOEXEPREFIX
     NOEXESUFFIX
     NAME exodus-test-field-utils
     COMM mpi serial
     NUM_MPI_PROCS 1
    )

  IF (NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")

#  TRIBITS_ADD_TEST(
#   testall
#   NOEXEPREFIX
#   NOEXESUFFIX
#   NAME exodus_unit_tests
#   COMM mpi serial
#   NUM_MPI_PROCS 1
#  )

  if ( TPL_Netcdf_Enables_Netcdf4 )
    TRIBITS_ADD_TEST(
     testall
     NOEXEPREFIX
     NOEXESUFFIX
     NAME exodus_unit_tests_nc4_env
     ARGS netcdf4
     COMM mpi serial
     NUM_MPI_PROCS 1
    )

    TRIBITS_ADD_TEST(
     testnc4
     NOEXEPREFIX
     NOEXESUFFIX
     NAME exodus_nc4_unit_tests
     COMM mpi serial
     NUM_MPI_PROCS 1
    )
  endif()

#  if ( TPL_Netcdf_Enables_Netcdf5 )
#    TRIBITS_ADD_TEST(
#      testall
#      NOEXEPREFIX
#      NOEXESUFFIX
#      NAME exodus_unit_tests_nc5_env
#      ARGS netcdf5
#      COMM mpi serial
#      NUM_MPI_PROCS 1
#    )
#  endif()
   endif()
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus unit tests will not be run due to failure finding ncdump executable.")
ENDIF()
