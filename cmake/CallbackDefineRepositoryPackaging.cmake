INCLUDE(TribitsGlobalMacros)

MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_REPOSITORY_DEFINE_PACKAGING() called for Trilinos!")
 
  # We need to make sure that these excludes only apply to Trilinos, not the global
  # project.
  IF (PROJECT_NAME STREQUAL Trilinos)
    SET(Trilinos_SOURCE_EXCLUDE_DIR "")
  ELSE()
    SET(Trilinos_SOURCE_EXCLUDE_DIR ${Trilinos_SOURCE_DIR})
  ENDIF()
  #PRINT_VAR(Trilinos_SOURCE_EXCLUDE_DIR)

  # Automatically update the version file for sierra
  IF (NOT MSVC)
    TRIBITS_REPOSITORY_CONFIGURE_VERSION_HEADER_FILE("Trilinos" "."
      ${Trilinos_SOURCE_DIR}/SIERRA/bjam/config_headers/${PROJECT_NAME}_version.h)
  ENDIF()

  SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    ${Trilinos_SOURCE_EXCLUDE_DIR}/.*[.]pyc$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/.*classicMakefile$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sparse_checkout.sh$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/tribits/common_tools/git/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/CMakeKitwareBacklog.txt$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/TODO$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/cmake/ctest/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/ITAPS/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/external/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/jpetra/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/cmmlib/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/configure.ac$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/configure$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/Makefile.am$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/Makefile.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/doc/[^b]
    ${Trilinos_SOURCE_EXCLUDE_DIR}/README_old$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sampleScripts/old_autotools/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/sampleScripts/git-profiles/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/SIERRA/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/coverage/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/harness/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/README$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/dependencies/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/packages/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/test/utilities/r.*
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/scripts/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/commonTools/release/
    ${Trilinos_SOURCE_EXCLUDE_DIR}/packages/common/DoxyfilePackageTemplate
    ${Trilinos_SOURCE_EXCLUDE_DIR}/stamp-h.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/configure.ac$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/aclocal.m4$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/configure$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/Makefile.am$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/Makefile.in$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/bootstrap$
    ${Trilinos_SOURCE_EXCLUDE_DIR}/config/
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
