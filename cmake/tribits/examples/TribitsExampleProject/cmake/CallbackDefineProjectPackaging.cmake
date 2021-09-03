macro(TRIBITS_PROJECT_DEFINE_PACKAGING)

    tribits_copy_installer_resource(TribitsExProj_README
      "${TribitsExProj_SOURCE_DIR}/README.md"
      "${TribitsExProj_BINARY_DIR}/README.md")
    tribits_copy_installer_resource(TribitsExProj_LICENSE
      "${TribitsExProj_SOURCE_DIR}/LICENSE"
      "${TribitsExProj_BINARY_DIR}/LICENSE.txt")

    set(CPACK_PACKAGE_DESCRIPTION "TribitsExampleProject just shows you how to use TriBITS.")
    set(CPACK_PACKAGE_FILE_NAME "tribitsexproj-setup-${TribitsExProj_VERSION}")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "TribitsExProj ${TribitsExProj_VERSION}")
    set(CPACK_PACKAGE_REGISTRY_KEY "TribitsExProj ${TribitsExProj_VERSION}")
    set(CPACK_PACKAGE_NAME "tribitsexproj")
    set(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
    set(CPACK_PACKAGE_VERSION "${TribitsExProj_VERSION}")
    set(CPACK_RESOURCE_FILE_README "${TribitsExProj_README}")
    set(CPACK_RESOURCE_FILE_LICENSE "${TribitsExProj_LICENSE}")
    set(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")
    set(CPACK_SOURCE_FILE_NAME "tribitsexproj-source-${TribitsExProj_VERSION}")
    set(CPACK_COMPONENTS_ALL ${TribitsExProj_PACKAGES} Unspecified)

endmacro()
