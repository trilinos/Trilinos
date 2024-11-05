include("${CMAKE_CURRENT_LIST_DIR}/../TribitsExProjCTestDriver.cmake")

set(COMM_TYPE SERIAL)
set(BUILD_TYPE DEBUG)
set(COMPILER_VERSION GCC)
set(BUILD_DIR_NAME ${COMM_TYPE}_${BUILD_TYPE})

set( EXTRA_CONFIGURE_OPTIONS
  "-DBUILD_SHARED_LIBS:BOOL=ON"
  "-DCMAKE_BUILD_TYPE=DEBUG"
  "-DCMAKE_C_COMPILER=gcc"
  "-DCMAKE_CXX_COMPILER=g++"
  "-DCMAKE_Fortran_COMPILER=gfortran"
  "-DTribitsExProj_ENABLE_Fortran=ON"
  "-DTribitsExProj_TRACE_ADD_TEST=ON"
  )

set_default_and_from_env(TribitsExProj_CMAKE_INSTALL_PREFIX "")
if (TribitsExProj_CMAKE_INSTALL_PREFIX)
  set(EXTRA_CONFIGURE_OPTIONS
    "${EXTRA_CONFIGURE_OPTIONS}"
    "-DCMAKE_INSTALL_PREFIX=${TribitsExProj_CMAKE_INSTALL_PREFIX}"
    )
endif()

set(CTEST_TEST_TYPE Continuous)

tribitsexproj_ctest_driver()
