  
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for s903186 using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME

  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET( CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  
  SET( CTEST_BUILD_FLAGS "-j8 -i" )

  SET_DEFAULT( CTEST_PARALLEL_LEVEL "8" )

  SET_DEFAULT( Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
  
  # Only turn on PyTrilinos for shared libraries
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES ${EXTRA_EXCLUDE_PACKAGES} TriKota Optika)
  
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DTrilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF"
    "-DBoost_INCLUDE_DIRS:FILEPATH=/Users/jmwille/install/boost_1_41_0"
    "-DTPL_ENABLE_Netcdf:BOOL=ON"
    "-DNetcdf_LIBRARY_DIRS=/Users/bmpersc/lib/netcdf-4.0/lib"
    "-DNetcdf_INCLUDE_DIRS=/Users/bmpersc/lib/netcdf-4.0/include"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE"
    "-DMesquite_ENABLE_TESTS:BOOL=ON"
    "-DCPPUNIT_LIBRARY:STRING=/Users/jmwille/install/cppunit-1.12.1/lib/libcppunit.a"
    "-DCPPUNIT_INCLUDES:STRING=/Users/jmwille/install/cppunit-1.12.1/include"
    )

  SET_DEFAULT(COMPILER_VERSION "GCC-4.6.1")
  
  IF (COMM_TYPE STREQUAL MPI)
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DTPL_ENABLE_MPI:BOOL=ON"
      "-DMPI_BASE_DIR:PATH=/Users/jmwille/install/gcc-4.6.1/openmpi-1.4.3"
      )
  
  ELSE()
  
    SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
      ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
      "-DCMAKE_CXX_COMPILER:FILEPATH=/Users/jmwille/install/gcc-4.6.1/bin/g++"
      "-DCMAKE_C_COMPILER:FILEPATH=/Users/jmwille/install/gcc-4.6.1/bin/gcc"
      "-DCMAKE_Fortran_COMPILER:FILEPATH=/Users/jmwille/install/gcc-4.6.1/bin/gfortran"
      )
  
  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
