INCLUDE(TribitsGlobalMacros)

MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_REPOSITORY_DEFINE_PACKAGING() called for Trilinos!")

  GET_FILENAME_COMPONENT(Trilinos_SOURCE_PATH ${Trilinos_SOURCE_DIR} PATH)

  # Automatically update the version file for sierra
  TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE("Trilinos" "."
    ${Trilinos_SOURCE_DIR}/SIERRA/bjam/config_headers/${PROJECT_NAME}_version.h)

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
    classicMakefile
    ".*.pyc"
    ${Trilinos_SOURCE_PATH}/cmake/tribits/common_tools/git
    ${Trilinos_SOURCE_PATH}/cmake/CMakeKitwareBacklog.txt
    ${Trilinos_SOURCE_PATH}/cmake/TODO
    ${Trilinos_SOURCE_PATH}/cmake/ctest
    ${Trilinos_SOURCE_PATH}/packages/ITAPS
    ${Trilinos_SOURCE_PATH}/packages/external
    ${Trilinos_SOURCE_PATH}/packages/jpetra
    ${Trilinos_SOURCE_PATH}/packages/cmmlib
    ${Trilinos_SOURCE_PATH}/packages/configure.ac
    ${Trilinos_SOURCE_PATH}/packages/configure
    ${Trilinos_SOURCE_PATH}/packages/Makefile.am
    ${Trilinos_SOURCE_PATH}/packages/Makefile.in
    ${Trilinos_SOURCE_PATH}/doc/[^b]
    ${Trilinos_SOURCE_PATH}/README_old
    ${Trilinos_SOURCE_PATH}/sampleScripts/old_autotools
    ${Trilinos_SOURCE_PATH}/sampleScripts/git-profiles
    ${Trilinos_SOURCE_PATH}/SIERRA
    ${Trilinos_SOURCE_PATH}/commonTools/test/coverage
    ${Trilinos_SOURCE_PATH}/commonTools/test/harness
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/README
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/dependencies
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/packages
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_PATH}/commonTools/scripts
    ${Trilinos_SOURCE_PATH}/commonTools/release
    ${Trilinos_SOURCE_PATH}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_PATH}/stamp-h.in
    ${Trilinos_SOURCE_PATH}/configure.ac
    ${Trilinos_SOURCE_PATH}/aclocal.m4
    ${Trilinos_SOURCE_PATH}/configure
    ${Trilinos_SOURCE_PATH}/Makefile.am
    ${Trilinos_SOURCE_PATH}/Makefile.in
    ${Trilinos_SOURCE_PATH}/bootstrap
    ${Trilinos_SOURCE_PATH}/config
  )
  
  #removing any packages not enabled from the tarball
  set(ENABLED_FLAG OFF)
  set(INCLUDE_EMPTY TRUE)
  TRIBITS_GET_ENABLED_LIST(${PROJECT_NAME}_PACKAGES ${PROJECT_NAME} ${ENABLED_FLAG} ${INCLUDE_EMPTY} 
    NON_ENABLED_PACKAGES NUM_NON_ENABLED)
  STRING(REPLACE " " ";" NON_ENABLED_PACKAGES "${NON_ENABLED_PACKAGES}")

  FOREACH(TRIBITS_PACKAGE ${NON_ENABLED_PACKAGES})
    #if the package is the TrilinosFramework we do not want to exclude it from the tarball
    #because that would exclude the cmake directory and the entire build system. So as a
    #special case we do not remove the TrilinosFramework from the tarball
    IF(NOT ${TRIBITS_PACKAGE} STREQUAL "TrilinosFramework")
      LIST(FIND ${PROJECT_NAME}_PACKAGES ${TRIBITS_PACKAGE} PACKAGE_IDX)
      LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
      
      #checking if we have a relative path to the package's files. Since the exclude is a
      #regular expression any "../" will be interpretted as <any char><any char>/ which
      #would never match the package's actual directory. There isn't a direct way in cmake
      #to convert a relative path into an absolute path with string operations so as a way
      #of making sure that we get the correct path of the package we use a find_path for the
      #CMakeLists.txt file for the package. Since the package has to have this file to work
      #correctly it should be guaranteed to be there.
      STRING(REGEX MATCH "[.][.]/" IS_RELATIVE_PATH ${PACKAGE_DIR})
      IF("${IS_RELATIVE_PATH}" STREQUAL "")
        SET(CPACK_SOURCE_IGNORE_FILES ${Trilinos_SOURCE_PATH}/${PACKAGE_DIR} ${CPACK_SOURCE_IGNORE_FILES})
      ELSE()
        FIND_PATH(ABSOLUTE_PATH CMakeLists.txt PATHS ${Trilinos_SOURCE_PATH}/${PACKAGE_DIR} NO_DEFAULT_PATH)
        IF("${ABSOLUTE_PATH}" STREQUAL "ABSOLUTE_PATH-NOTFOUND")
          MESSAGE(AUTHOR_WARNING "Relative path found for disabled package ${TRIBITS_PACKAGE} but package was missing a CMakeLists.txt file. This disabled package will likely not be excluded from a source release")
        ENDIF()
        SET(CPACK_SOURCE_IGNORE_FILES ${ABSOLUTE_PATH} ${CPACK_SOURCE_IGNORE_FILES})
      ENDIF()
    ENDIF()
  ENDFOREACH()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("Exclude files when building source packages")
    FOREACH(item ${CPACK_SOURCE_IGNORE_FILES})
      MESSAGE(${item})
    ENDFOREACH()
  ENDIF()

  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the pertinant file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  MACRO(COPY_INSTALLER_RESOURCE _varname _source _destination)
    SET("${_varname}" "${_destination}")
    IF (EXISTS "${_destination}")
      FILE(REMOVE_RECURSE "${_destination}")
    ENDIF ()
    CONFIGURE_FILE(
      "${_source}" 
      "${_destination}" 
      COPYONLY)
  ENDMACRO()
  COPY_INSTALLER_RESOURCE(Trilinos_README
    "${Trilinos_SOURCE_DIR}/README"
    "${Trilinos_BINARY_DIR}/README.txt")
  COPY_INSTALLER_RESOURCE(Trilinos_LICENSE
    "${Trilinos_SOURCE_DIR}/LICENSE"
    "${Trilinos_BINARY_DIR}/LICENSE.txt")

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_RESOURCE_FILE_README "${Trilinos_README}")
  SET(CPACK_RESOURCE_FILE_LICENSE "${Trilinos_LICENSE}")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES} Unspecified)
  
  set(ENABLED_FLAG ON)
  set(INCLUDE_EMPTY FALSE)
  TRIBITS_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ${ENABLED_FLAG}
    ${INCLUDE_EMPTY} ENABLED_PACKAGES NUM_ENABLED)
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
    #Resetting the name to avoid overwriting registery keys when installing
    SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-${Trilinos_VERSION}")
    IF (TPL_ENABLE_MPI)
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-mpi")
    ELSE ()
      SET(CPACK_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-serial")
    ENDIF()
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH OFF)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()
