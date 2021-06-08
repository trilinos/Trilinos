# Copyright(C) 1999-2021 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

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
    testrdd
    testwtdi64
    testrddi64
  )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/testall.in
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall
    @ONLY
  )
  # now copy the temporary into the final destination, setting the permissions
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  # NOTE: The "LINKER_LANGUAGE CXX" is used to workaround an issue with
  #       the trilinos cmake adding the pamgen and pamgen_extras libraries as
  #       a dependency.  These are C++ libraries and result in unresolved dependencies
  #       when built with shareable libraries on some platforms....

  FOREACH ( testName ${forbind_OLDTESTS} )
    TRIBITS_ADD_EXECUTABLE( f32_${testName} NOEXEPREFIX NOEXESUFFIX LINKER_LANGUAGE Fortran SOURCES ${testName}.f)
  ENDFOREACH ( testName )

  TRIBITS_ADD_TEST(
   testall
   NOEXEPREFIX
   NOEXESUFFIX
   NAME exodus_unit_tests_for32
   COMM mpi serial
   NUM_MPI_PROCS 1
  )

  if ( TPL_Netcdf_Enables_Netcdf4 )
     CONFIGURE_FILE(
       ${CMAKE_CURRENT_SOURCE_DIR}/testall_i64.in
       ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall_i64
       @ONLY
     )
     # now copy the temporary into the final destination, setting the permissions
     file(COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/testall_i64
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

    TRIBITS_ADD_TEST(
      testall
      NOEXEPREFIX
      NOEXESUFFIX
      NAME exodus_unit_tests_for32_64
      ARGS netcdf4
      COMM mpi serial
      NUM_MPI_PROCS 1
    )

    TRIBITS_ADD_TEST(
      testall_i64
      NOEXEPREFIX
      NOEXESUFFIX
      NAME exodus_nc4_unit_tests_for32
      COMM mpi serial
      NUM_MPI_PROCS 1
    )
  endif()
ELSE()
  MESSAGE(STATUS
      "WARNING: ${PACKAGE_NAME}: Exodus unit tests will not be run due to failure finding ncdump executable.")
ENDIF()
