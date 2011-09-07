

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
    ${Trilinos_SOURCE_DIR}/packages/seacas
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

  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the README file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  SET(Trilinos_LICENSE_AND_README "${Trilinos_BINARY_DIR}/README.txt")
  FILE(COPY "${Trilinos_SOURCE_DIR}/README" DESTINATION ${Trilinos_LICENSE_AND_README})

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_LICENSE_AND_README}")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_LICENSE_AND_README}")
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
