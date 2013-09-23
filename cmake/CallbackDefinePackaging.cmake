INCLUDE(TribitsGlobalMacros)

MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_REPOSITORY_DEFINE_PACKAGING() called for Trilinos!")

  GET_FILENAME_COMPONENT(Trilinos_SOURCE_PATH ${Trilinos_SOURCE_DIR} PATH)

  # Automatically update the version file for sierra
  IF (NOT MSVC)
    TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE("Trilinos" "."
      ${Trilinos_SOURCE_DIR}/SIERRA/bjam/config_headers/${PROJECT_NAME}_version.h)
  ENDIF()

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
    classicMakefile
    ${Trilinos_SOURCE_PATH}/sparse_checkout.sh
    ${Trilinos_SOURCE_PATH}/cmake/tribits/common_tools/git/
    ${Trilinos_SOURCE_PATH}/cmake/CMakeKitwareBacklog.txt
    ${Trilinos_SOURCE_PATH}/cmake/TODO
    ${Trilinos_SOURCE_PATH}/cmake/ctest/
    ${Trilinos_SOURCE_PATH}/packages/ITAPS/
    ${Trilinos_SOURCE_PATH}/packages/external/
    ${Trilinos_SOURCE_PATH}/packages/jpetra/
    ${Trilinos_SOURCE_PATH}/packages/cmmlib/
    ${Trilinos_SOURCE_PATH}/packages/configure.ac
    ${Trilinos_SOURCE_PATH}/packages/configure
    ${Trilinos_SOURCE_PATH}/packages/Makefile.am
    ${Trilinos_SOURCE_PATH}/packages/Makefile.in
    ${Trilinos_SOURCE_PATH}/doc/[^b]
    ${Trilinos_SOURCE_PATH}/README_old
    ${Trilinos_SOURCE_PATH}/sampleScripts/old_autotools/
    ${Trilinos_SOURCE_PATH}/sampleScripts/git-profiles/
    ${Trilinos_SOURCE_PATH}/SIERRA/
    ${Trilinos_SOURCE_PATH}/commonTools/test/coverage/
    ${Trilinos_SOURCE_PATH}/commonTools/test/harness/
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/README
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/dependencies/
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/packages/
    ${Trilinos_SOURCE_PATH}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_PATH}/commonTools/scripts/
    ${Trilinos_SOURCE_PATH}/commonTools/release/
    ${Trilinos_SOURCE_PATH}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_PATH}/stamp-h.in
    ${Trilinos_SOURCE_PATH}/configure.ac
    ${Trilinos_SOURCE_PATH}/aclocal.m4
    ${Trilinos_SOURCE_PATH}/configure
    ${Trilinos_SOURCE_PATH}/Makefile.am
    ${Trilinos_SOURCE_PATH}/Makefile.in
    ${Trilinos_SOURCE_PATH}/bootstrap
    ${Trilinos_SOURCE_PATH}/config/
    )

  SET(TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE ${TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE}
    TrilinosFramework)

  IF (PROJECT_NAME STREQUAL Trilinos)
  
    # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
    # .txt .rtf .html. Copy the pertinant file to the binary directory with
    # a .txt extension. This is only the case with the PackageMaker 
    # generator, but it doesn't hurt to do it for other generators as
    # well.
    TRIBITS_COPY_INSTALLER_RESOURCE(Trilinos_README
      "${Trilinos_SOURCE_DIR}/README"
      "${Trilinos_BINARY_DIR}/README.txt")
    TRIBITS_COPY_INSTALLER_RESOURCE(Trilinos_LICENSE
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
#    SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
    SET(CPACK_SOURCE_GENERATOR "TGZ")
    SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
    SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES} Unspecified)

  ENDIF()  
  
ENDMACRO()
