INCLUDE(Join)
INCLUDE(TribitsFilepathHelpers)


FUNCTION(TRIBITS_CONFIGURE_SET_INSTALLED_GROUP_AND_PERMS_FILE  TARGET_FILE)

  SET(PROJECT_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR
    "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}")

  TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR(
    "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}"
    "${CMAKE_INSTALL_PREFIX}"
    PROJECT_SUBDIR_PATHS_ARRAY
    )

  SET(PROJECT_MAKE_INSTALL_GROUP "${${PROJECT_NAME}_MAKE_INSTALL_GROUP}")

  SET(group_perms "")
  IF (${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE)
    SET(group_perms "g+rwX")
  ELSEIF (${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE)
    SET(group_perms "g+rX")
  ENDIF()

  SET(other_perms "")
  IF (${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE)
    SET(other_perms "o+rX")
  ENDIF()

  JOIN(PROJECT_MAKE_INSTALL_PERMS_CHANGE "," FALSE
    ${group_perms} ${other_perms} )

  SET(tribits_install_src
    "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")
  CONFIGURE_FILE(
    "${tribits_install_src}/set_installed_group_and_permissions.cmake.in"
    "${TARGET_FILE}" @ONLY )

ENDFUNCTION()


FUNCTION(TRIBITS_ADD_INSTALL_GROUP_AND_PERMS_FIXUPS)

  IF (NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")

    SET(set_installed_group_and_permissions_file
      "${PROJECT_BINARY_DIR}/set_installed_group_and_permissions.cmake")

    TRIBITS_CONFIGURE_SET_INSTALLED_GROUP_AND_PERMS_FILE(
      "${set_installed_group_and_permissions_file}" )

    # Fix up install for default 'install' command
    INSTALL(SCRIPT "${set_installed_group_and_permissions_file}")

  ENDIF()

ENDFUNCTION()