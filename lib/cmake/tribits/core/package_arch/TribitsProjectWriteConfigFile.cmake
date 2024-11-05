# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
#
# Module TribitsInternalPackageWriteConfigFile.cmake
#
# This module contains code for generating <Project>Config.cmake files for
# TriBITS projects.
#
################################################################################


include(TribitsInternalPackageWriteConfigFile)


# Write the outer TriBITS <Project>Config.cmake file
#
# If ${PROJECT_NAME}_VERSION is not set or is '' on input, then it will be set
# to 0.0.0 in order to create the ${PROJECT_NAME}ConfigVersion.cmake file.
#
# ToDo: Finish documentation!
#
function(tribits_write_project_client_export_files)

  set(EXPORT_FILE_VAR_PREFIX ${PROJECT_NAME})

  # Reversing the package list so that libraries will be produced in order of
  # most dependent to least dependent.
  set(PACKAGE_LIST ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
  if (PACKAGE_LIST)
    list(REVERSE PACKAGE_LIST)
  endif()

  # Loop over all packages to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  set(FULL_PACKAGE_SET "")
  set(FULL_LIBRARY_SET "")
  foreach(TRIBITS_PACKAGE ${PACKAGE_LIST})
    if(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      list(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      list(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
    endif()
  endforeach()

  set(${PROJECT_NAME}_CONFIG_LIBRARIES ${FULL_LIBRARY_SET})

  # Reversing the tpl list so that the list of tpls will be produced in
  # order of most dependent to least dependent.
  if (${PROJECT_NAME}_DEFINED_TPLS)
    set(TPL_LIST ${${PROJECT_NAME}_DEFINED_TPLS})
    list(REVERSE TPL_LIST)
  endif()

  # Loop over all TPLs to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  set(FULL_TPL_SET "")
  set(FULL_TPL_LIBRARY_SET "")
  foreach(TPL ${TPL_LIST})
    if(TPL_ENABLE_${TPL})
      list(APPEND FULL_TPL_SET ${TPL})
      list(APPEND FULL_TPL_LIBRARY_SET ${TPL_${TPL}_LIBRARIES})
    endif()
  endforeach()

  set(${PROJECT_NAME}_CONFIG_TPL_LIBRARIES ${FULL_TPL_LIBRARY_SET})

  #
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree
  # and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  set(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    string(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND
     "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
    set(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  endif()

  # Custom code in configuration file.
  set(PROJECT_CONFIG_CODE "")

  # Appending the logic to include each package's config file.
  set(LOAD_CODE "# Load configurations from enabled packages")
  foreach(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    set(LOAD_CODE "${LOAD_CODE}
include(\"${${TRIBITS_PACKAGE}_BINARY_DIR}/${TRIBITS_PACKAGE}Config.cmake\")")
  endforeach()
  set(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  tribits_set_compiler_vars_for_config_file(INSTALL_DIR)

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    # In TribitsProjectConfigTemplate.cmake.in, we would like to preserve
    # ${}-variables after the conversion to TribitsProjectConfigTemplate.cmake.
    # To this end, one typically uses the @-syntax for variables. That doesn't
    # support nested variables, however.  Use ${PDOLLAR} as a workaround, cf.
    # <http://www.cmake.org/pipermail/cmake/2013-April/054341.html>.
    set(PDOLLAR "$")
    set(tribitsInstallationDir
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")
    set(TRIBITS_PROJECT_INSTALL_INCLUDE_DIR "")
    configure_file(
      "${tribitsInstallationDir}/TribitsProjectConfigTemplate.cmake.in"
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" )
  endif()

  ######
  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree. It will be placed in the base
  # install directory for ${PROJECT_NAME} when installed.
  ######

  # Set the include and library directories relative to the location
  # at which the ${PROJECT_NAME}Config.cmake file is going to be
  # installed. Note the variable reference below is escaped so it
  # won't be replaced until a client project attempts to locate
  # directories using the installed config file. This is to deal with
  # installers that allow relocation of the install tree at *install*
  # time.
  # The export files are typically installed in
  #     <install dir>/<lib path>/cmake/<package name>/.
  # The relative path to the installation dir is hence k*(../) + ../../, where
  # k is the number of components in <lib path>. Extract those here.
  # This doesn't work if ${${PROJECT_NAME}_INSTALL_LIB_DIR} contains "./" or
  # "../" components, but really, it never did. All of this should actually be
  # handled by CMake's configure_package_config_file().
  string(REPLACE "/" ";" PATH_LIST ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  set(RELATIVE_PATH "../..")
  foreach(PATH ${PATH_LIST})
    set(RELATIVE_PATH "${RELATIVE_PATH}/..")
  endforeach()

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    set(SHARED_LIB_RPATH_COMMAND
       "${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}"
      )
  endif()

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)

    tribits_set_compiler_vars_for_config_file(INSTALL_DIR)

    # Custom code in configuration file.
    set(PROJECT_CONFIG_CODE "")

    set(PDOLLAR "$")  # Hack used in configure file below

    if (IS_ABSOLUTE "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
      set(TRIBITS_PROJECT_INSTALL_INCLUDE_DIR "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
    else()
      set(TRIBITS_PROJECT_INSTALL_INCLUDE_DIR
        "${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
    endif()

    configure_file(
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in"
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake"
      )

    install(
      FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake"
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      RENAME ${PROJECT_NAME}Config.cmake
      )
  endif()

  #
  # Configure the version file for ${PROJECT_NAME}
  #
  include(CMakePackageConfigHelpers)
  if ("${${PROJECT_NAME}_VERSION}"  STREQUAL  "")
    set(${PROJECT_NAME}_VERSION  0.0.0)
  endif()
  write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${${PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion
    )
  install(
    FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
    )

endfunction()
