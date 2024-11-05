include(Join)
include(TribitsFilepathHelpers)
include(AppendStringVarWithSep)


function(tribits_raise_install_perms_mods_not_supported_on_windows_error)

  set(INSTALL_PERMS_SET "")
  tribits_append_install_perms_var_not_supported(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE)
  tribits_append_install_perms_var_not_supported(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE)
  tribits_append_install_perms_var_not_supported(
    ${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE)
  tribits_append_install_perms_var_not_supported(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP)

  message(FATAL_ERROR
    "ERROR: The options:\n"
    "${INSTALL_PERMS_SET}"
    "are not supported on Windows!\n"
    "Please remove these options and configure from scratch!"
    )

endfunction()


# Reads and writes var INSTALL_PERMS_SET in above function
macro(tribits_append_install_perms_var_not_supported  VAR_NAME)
  if (NOT "${${VAR_NAME}}" STREQUAL "")
    set(INSTALL_PERMS_SET  "${INSTALL_PERMS_SET}    ${VAR_NAME}='${${VAR_NAME}}'\n")
  endif()
endmacro()


function(tribits_determine_if_setup_for_group_and_perms_modifications
  SETUP_FOR_GROUP_AND_PERMS_MODIFICATIONS_OUT
  )

  if(
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE OR
    ${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE OR
    ${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE OR
    (NOT "${${PROJECT_NAME}_MAKE_INSTALL_GROUP}" STREQUAL "")
    )
    set(setupForGroupAndPermsModifications TRUE)
  else()
    set(setupForGroupAndPermsModifications FALSE)
  endif()

  if (setupForGroupAndPermsModifications AND
    ${PROJECT_NAME}_HOSTTYPE STREQUAL "Windows"
    )
    tribits_raise_install_perms_mods_not_supported_on_windows_error()
  endif()

  set(${SETUP_FOR_GROUP_AND_PERMS_MODIFICATIONS_OUT}
    ${setupForGroupAndPermsModifications} PARENT_SCOPE)

endfunction()


function(tribits_configure_set_installed_group_and_perms_file  TARGET_FILE)

  set(PROJECT_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR
    "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}")

  tribits_get_dir_array_below_base_dir(
    "${${PROJECT_NAME}_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR}"
    "${CMAKE_INSTALL_PREFIX}"
    PROJECT_SUBDIR_PATHS_ARRAY
    )

  set(PROJECT_MAKE_INSTALL_GROUP "${${PROJECT_NAME}_MAKE_INSTALL_GROUP}")

  set(group_perms "")
  if (${PROJECT_NAME}_MAKE_INSTALL_GROUP_WRITABLE)
    set(group_perms "g+rwX")
  elseif (${PROJECT_NAME}_MAKE_INSTALL_GROUP_READABLE
    OR ${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE
    )
    set(group_perms "g+rX")
  endif()

  set(other_perms "")
  if (${PROJECT_NAME}_MAKE_INSTALL_WORLD_READABLE)
    set(other_perms "o+rX")
  endif()

  join(PROJECT_MAKE_INSTALL_PERMS_CHANGE "," FALSE
    ${group_perms} ${other_perms} )

  set(tribits_install_src
    "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")
  configure_file(
    "${tribits_install_src}/set_installed_group_and_permissions.cmake.in"
    "${TARGET_FILE}" @ONLY )

endfunction()


function(tribits_add_install_group_and_perms_fixups)

  tribits_determine_if_setup_for_group_and_perms_modifications(
    setupForGroupAndPermsModifications)

  if (setupForGroupAndPermsModifications)

    set(set_installed_group_and_permissions_file
      "${PROJECT_BINARY_DIR}/set_installed_group_and_permissions.cmake")

    tribits_configure_set_installed_group_and_perms_file(
      "${set_installed_group_and_permissions_file}" )

    # Fix up install for default 'install' command
    install(SCRIPT "${set_installed_group_and_permissions_file}")

  endif()

endfunction()