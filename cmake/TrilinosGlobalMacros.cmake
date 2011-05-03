

#
# This file contains global-level macros that are specific to Trilinos
#



#
# Macro that defines Trilinos testing support
#

MACRO(TRILINOS_SETUP_TESTING_SUPPORT)

  CONFIGURE_FILE(
    ${Trilinos_SOURCE_DIR}/cmake/ctest/CTestCustom.ctest.in
    ${Trilinos_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()


#
# Macro that drives a Trilinos 'dashboard' target
#

MACRO(TRILINOS_ADD_DASHBOARD_TARGET)

  IF (NOT (WIN32 AND NOT CYGWIN))

    ADVANCED_SET(Trilinos_DASHBOARD_CTEST_ARGS "" CACHE STRING
      "Extra arguments to pass to CTest when calling 'ctest -S' to run the 'dashboard' make target." )
  
    # H.1) Enable all packages that are enabled and have tests enabled
  
    SET(Trilinos_ENABLED_PACKAGES_LIST)
    SET(Trilinos_ENABLED_PACKAGES_CMAKE_ARG_LIST)
    FOREACH(PACKAGE ${Trilinos_PACKAGES})
      IF (Trilinos_ENABLE_${PACKAGE} AND ${PACKAGE}_ENABLE_TESTS)
        IF (Trilinos_ENABLED_PACKAGES_LIST)
          SET(Trilinos_ENABLED_PACKAGES_LIST
            "${Trilinos_ENABLED_PACKAGES_LIST}\;${PACKAGE}") 
        ELSE()
          SET(Trilinos_ENABLED_PACKAGES_LIST "${PACKAGE}") 
        ENDIF()
        SET(Trilinos_ENABLED_PACKAGES_CMAKE_ARG_LIST
          ${Trilinos_ENABLED_PACKAGES_CMAKE_ARG_LIST} -DTrilinos_ENABLE_${PACKAGE}=ON)
      ENDIF()
    ENDFOREACH()
    #PRINT_VAR(Trilinos_ENABLED_PACKAGES_LIST)
    
    SET(EXPR_CMND_ARGS)
    IF (Trilinos_ENABLE_COVERAGE_TESTING)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DO_COVERAGE_TESTING=TRUE")
    ENDIF()

    #PRINT_VAR(${PROJECT_NAME}_EXTRA_REPOSITORIES)
    JOIN(${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED "," FALSE
      ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRA_REPOSITORIES=${${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED})
  
    # H.2) Add the custom target to enable all the packages with tests enabled
    
    ADD_CUSTOM_TARGET(dashboard
  
      VERBATIM
    
      # WARNING: The echoed command and the actual commands are duplicated!  You have to reproduce them!
  
      COMMAND echo
      COMMAND echo "***************************************************"
      COMMAND echo "*** Running incremental experimental dashboard ***" 
      COMMAND echo "***************************************************"
      COMMAND echo
      COMMAND echo Trilinos_ENABLED_PACKAGES_LIST=${Trilinos_ENABLED_PACKAGES_LIST}
      COMMAND echo
  
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** A) Clean out the list of packages"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} -DTrilinos_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF ${TRILINOS_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} -DTrilinos_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF ${TRILINOS_HOME_DIR}
  
      # NOTE: Above, if Trilinos_ENABLE_ALL_PACKAGES was set in CMakeCache.txt, then setting
      # -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF will turn it off in the cache.  Note that it will
      # never be turned on again which means that the list of packages will be set explicitly below.
    
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** B) Run the dashboard command setting the list of packages"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: env ${EXPR_CMND_ARGS}
        Trilinos_PACKAGES=${Trilinos_ENABLED_PACKAGES_LIST}
        ${CMAKE_CTEST_COMMAND} ${Trilinos_DASHBOARD_CTEST_ARGS} -S
          ${TRILINOS_HOME_DIR}/cmake/ctest/experimental_build_test.cmake
      COMMAND echo
      COMMAND env ${EXPR_CMND_ARGS}
        Trilinos_PACKAGES=${Trilinos_ENABLED_PACKAGES_LIST}
        ${CMAKE_CTEST_COMMAND} ${Trilinos_DASHBOARD_CTEST_ARGS} -S
          ${TRILINOS_HOME_DIR}/cmake/ctest/experimental_build_test.cmake || echo
  
      # 2009/07/05: rabartl: Above, I added the ending '|| echo' to always make
      # the command pass so that 'make will not stop and avoid this last command
      # to set back the enabled packages.
  
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** C) Clean out the list of packages again to clean the cache file"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} -DTrilinos_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF ${TRILINOS_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} -DTrilinos_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF ${TRILINOS_HOME_DIR}
    
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** D) Reconfigure with the original package list"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} ${Trilinos_ENABLED_PACKAGES_CMAKE_ARG_LIST}
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON ${TRILINOS_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} ${Trilinos_ENABLED_PACKAGES_CMAKE_ARG_LIST}
        -DTrilinos_ALLOW_NO_PACKAGES:BOOL=ON ${TRILINOS_HOME_DIR}
  
      COMMAND echo
      COMMAND echo "See the results at http://trilinos.sandia.gov/cdash/index.php?project=Trilinos&display=project\#Experimental"
      COMMAND echo
   
      )
  
  ENDIF()

ENDMACRO()


#
# Macro that defines Trilinos packaging options:
#

MACRO(TRILINOS_DEFINE_PACKAGING)

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
    classicMakefile
    ".*.pyc"
    ${Trilinos_SOURCE_DIR}/cmake/CMakeKitwareBacklog.txt
    ${Trilinos_SOURCE_DIR}/cmake/TODO
    ${Trilinos_SOURCE_DIR}/packages/ITAPS
    ${Trilinos_SOURCE_DIR}/packages/aristos
    ${Trilinos_SOURCE_DIR}/packages/claps
    ${Trilinos_SOURCE_DIR}/packages/external
    ${Trilinos_SOURCE_DIR}/packages/jpetra
    ${Trilinos_SOURCE_DIR}/packages/new_package
    ${Trilinos_SOURCE_DIR}/packages/rbgen
    ${Trilinos_SOURCE_DIR}/packages/WebTrilinos
    ${Trilinos_SOURCE_DIR}/packages/cmmlib
    ${Trilinos_SOURCE_DIR}/packages/Trios
    ${Trilinos_SOURCE_DIR}/packages/meros
    ${Trilinos_SOURCE_DIR}/packages/phdmesh
    ${Trilinos_SOURCE_DIR}/demos/FEApp
    ${Trilinos_SOURCE_DIR}/packages/configure.ac
    ${Trilinos_SOURCE_DIR}/packages/configure
    ${Trilinos_SOURCE_DIR}/packages/Makefile.am
    ${Trilinos_SOURCE_DIR}/packages/Makefile.in
    ${Trilinos_SOURCE_DIR}/doc/[^b]
    ${Trilinos_SOURCE_DIR}/README_old
    ${Trilinos_SOURCE_DIR}/sampleScripts/old_autotools
    ${Trilinos_SOURCE_DIR}/sampleScripts/git-profiles
    ${Trilinos_SOURCE_DIR}/SIERRA/
    ${Trilinos_SOURCE_DIR}/commonTools/test/coverage
    ${Trilinos_SOURCE_DIR}/commonTools/test/harness
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/README
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/dependencies
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/packages
    ${Trilinos_SOURCE_DIR}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_DIR}/commonTools/scripts
    ${Trilinos_SOURCE_DIR}/commonTools/release
    ${Trilinos_SOURCE_DIR}/commonTools/git
    ${Trilinos_SOURCE_DIR}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_DIR}/stamp-h.in
    ${Trilinos_SOURCE_DIR}/configure.ac
    ${Trilinos_SOURCE_DIR}/aclocal.m4
    ${Trilinos_SOURCE_DIR}/configure
    ${Trilinos_SOURCE_DIR}/Makefile.am
    ${Trilinos_SOURCE_DIR}/Makefile.in
    ${Trilinos_SOURCE_DIR}/bootstrap
    ${Trilinos_SOURCE_DIR}/config
  )

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()
  

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_SOURCE_DIR}/README")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES})
  
  PACKAGE_ARCH_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ON
    FALSE ENABLED_PACKAGES NUM_ENABLED)
  string(REPLACE " " ";" ENABLED_PACKAGES "${ENABLED_PACKAGES}")
  
  #message("ENABLED PACKAGES: ${ENABLED_PACKAGES} ${NUM_ENABLED}")
  FOREACH(PKG ${ENABLED_PACKAGES})
    IF(NOT "${${PKG}_LIB_REQUIRED_DEP_PACKAGES}" STREQUAL "")
        string(TOUPPER ${PKG} UPPER_PKG)
        #message("${UPPER_PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
        SET(CPACK_COMPONENT_${UPPER_PKG}_DEPENDS ${${PKG}_LIB_REQUIRED_DEP_PACKAGES})
    ENDIF()
    #message("${PKG} depends on : ${${PKG}_LIB_REQUIRED_DEP_PACKAGES}")
  ENDFOREACH()

  
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
