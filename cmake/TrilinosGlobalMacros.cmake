

#
# This file contains global-level macros that are specific to Trilinos
#



#
# Macro that defines Trilinos testing support
#

MACRO(TRILINOS_SETUP_TESTING_SUPPORT)
  
  IF (WIN32 AND NOT CYGWIN)
    SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT OFF)
  ELSE()
    SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT ON)
  ENDIF()
  
  ADVANCED_OPTION(Trilinos_ENABLE_NATIVE_TEST_HARNESS
    "Enable the native Trilinos perl-based test harness."
    ${Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT} )
  
  IF (Trilinos_ENABLE_NATIVE_TEST_HARNESS)
  
    ADD_CUSTOM_TARGET(
      runtests-serial
       ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
      --trilinos-dir=${TRILINOS_HOME_DIR}
      --comm=serial
      --build-dir=${TRILINOS_BUILD_DIR}
      --category=${TRILINOS_TEST_CATEGORY}
      --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
      )
  
    IF (TPL_ENABLE_MPI)
    
      ADD_CUSTOM_TARGET(
        runtests-mpi
         ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
        --trilinos-dir=${TRILINOS_HOME_DIR}
        --comm=mpi
        --mpi-go="${TRILINOS_MPI_GO}"
        --max-proc=${MPIEXEC_MAX_NUMPROCS}
        --build-dir=${TRILINOS_BUILD_DIR}
        --category=${TRILINOS_TEST_CATEGORY}
        --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
        )
  
    ENDIF()
  
  ENDIF()
  
  IF (WIN32)
    SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT OFF)
  ELSE()
    SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT ON)
  ENDIF()
  
  # 2008/10/17: rabartl: Above, I can not turn these tests on by default
  # with cygwin because the custom script target is not working for some
  # reason.
  
  ADVANCED_OPTION(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS
    "Enable dependency unit tests."
    ${Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT}
    )

  # Add the directory for the unit tests
  ADD_SUBDIRECTORY(cmake)

  CONFIGURE_FILE(
    ${Trilinos_SOURCE_DIR}/cmake/ctest/CTestCustom.ctest.in
    ${Trilinos_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()


#
# Macro that defines Trilinos packaging options:
#

MACRO(TRILINOS_DEFINE_PACKAGING)

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH ON)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
