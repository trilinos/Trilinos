MACRO(TRIBITS_PROJECT_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_PROJECT_DEFINE_PACKAGING() called for Trilinos!")
  
  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the pertinant file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  TRIBITS_COPY_INSTALLER_RESOURCE(Trilinos_README
    "${Trilinos_SOURCE_DIR}/README.md"
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
    SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")
    SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
    SET(CPACK_COMPONENTS_ALL ${Trilinos_PACKAGES} Unspecified)
  
ENDMACRO()
