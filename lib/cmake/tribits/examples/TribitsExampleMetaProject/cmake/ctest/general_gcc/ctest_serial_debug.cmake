include("${CMAKE_CURRENT_LIST_DIR}/../TribitsExMetaProjCTestDriver.cmake")

set(COMM_TYPE SERIAL)
set(BUILD_TYPE DEBUG)
set(COMPILER_VERSION GCC)
set(BUILD_DIR_NAME ${COMM_TYPE}_${BUILD_TYPE})

set_default( CTEST_BUILD_FLAGS "-j1 -i" )
set_default( CTEST_PARALLEL_LEVEL "1" )

set( EXTRA_CONFIGURE_OPTIONS
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DCMAKE_BUILD_TYPE=DEBUG"
  "-DCMAKE_C_COMPILER=gcc"
  "-DCMAKE_CXX_COMPILER=g++"
  "-DCMAKE_Fortran_COMPILER=gfortran"
  "-DTribitsExMetaProj_ENABLE_Fortran=ON"
  "-DTribitsExMetaProj_TRACE_ADD_TEST=ON"
  )

set_default_and_from_env(TribitsExMetaProj_CMAKE_INSTALL_PREFIX "")
if (TribitsExMetaProj_CMAKE_INSTALL_PREFIX)
  set(EXTRA_CONFIGURE_OPTIONS
    "${EXTRA_CONFIGURE_OPTIONS}"
    "-DCMAKE_INSTALL_PREFIX=${TribitsExMetaProj_CMAKE_INSTALL_PREFIX}"
    )
endif()

set(CTEST_TEST_TYPE Continuous)

tribitsexmetaproj_ctest_driver()
