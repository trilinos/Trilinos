
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for IBM AIX machine p90n03
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)
  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  SET( CTEST_DASHBOARD_ROOT "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )

  SET(ENV{OBJECT_MODE} 64)
  SET(COMPILER_VERSION XL-11.01)
  SET(CTEST_NOTES_FILES "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  SET(CTEST_BUILD_FLAGS "-j8 -i" )
  SET(Trilinos_REPOSITORY_LOCATION $ENV{HOME}/bench/gitroot/Trilinos.git)
  SET(Trilinos_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE " ")
  SET(Trilinos_WARNINGS_AS_ERRORS_FLAGS " ")
  #SET(Trilinos_TRACK "Experimental")
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES PyTrilinos TriKota Optika)
  SET( EXTRA_SYSTEM_CONFIGURE_OPTIONS
    "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"
    "-DCMAKE_C_COMPILER:FILEPATH=/usr/vacpp/bin/xlc_r"
    "-DCMAKE_CXX_COMPILER:FILEPATH=/usr/vacpp/bin/xlC_r"
    "-DCMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/xlf_r"
    "-DCMAKE_CXX_FLAGS:STRING=-qrtti=all -qstaticinline"
    "-DCMAKE_Fortran_FLAGS:STRING=-qarch=pwr7 -qtune=pwr7"
    "-DBUILD_SHARED_LIBS=ON"
    "-DTrilinos_ENABLE_TriKota:BOOL=OFF"
    "-DTPL_BLAS_LIBRARIES=$ENV{LAPACK_DIR}/libblas.so"
    "-DTPL_LAPACK_LIBRARIES=$ENV{LAPACK_DIR}/liblapack.so"
    "-DTPL_Pthread_LIBRARIES=-lpthread"
    )

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
