SET(tribits_install_src
  "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")

IF((${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES
    OR ${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
  AND NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
  )

  INCLUDE(TribitsWriteClientExportFiles)

  TRIBITS_WRITE_PROJECT_CLIENT_EXPORT_FILES()

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    # TEMPORARY: Install a compatibility copy of ${PROJECT_NAME}Config.cmake
    # where was previously installed to warn and load the new file.
    SET(COMPATIBILITY_CONFIG_INCLUDE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
    CONFIGURE_FILE(
      "${tribits_install_src}/TribitsConfigInclude.cmake.in"
      ${COMPATIBILITY_CONFIG_INCLUDE}
      @ONLY
      )
    INSTALL(
      FILES ${COMPATIBILITY_CONFIG_INCLUDE}
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      )
  ENDIF()

ENDIF()

# NOTE: This install logic is done in this subdir CMakeList.txt file so that
# the install_package_by_package install target can invoke it for this subdir
# even if any of the package installs fail!  See the script
# 'cmake_pbp_install.cmake.in' for more details.
