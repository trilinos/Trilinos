set(tribits_install_src
  "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")

if ( ${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
    AND (NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    AND (NOT ${PROJECT_NAME}_SKIP_INSTALL_PROJECT_CMAKE_CONFIG_FILES)
  )

  include(TribitsProjectWriteConfigFile)

  tribits_write_project_client_export_files()

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    # TEMPORARY: Install a compatibility copy of ${PROJECT_NAME}Config.cmake
    # where was previously installed to warn and load the new file.
    set(COMPATIBILITY_CONFIG_INCLUDE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
    configure_file(
      "${tribits_install_src}/TribitsConfigInclude.cmake.in"
      ${COMPATIBILITY_CONFIG_INCLUDE}
      @ONLY
      )
    install(
      FILES ${COMPATIBILITY_CONFIG_INCLUDE}
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      )
  endif()

endif()

# NOTE: This install logic is done in this subdir CMakeList.txt file so that
# the install_package_by_package install target can invoke it for this subdir
# even if any of the package installs fail!  See the script
# 'cmake_pbp_install.cmake.in' for more details.
