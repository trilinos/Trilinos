# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

include(TribitsGeneralMacros)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake!
###

#  This function will take a list and turn it into a space separated string
#  adding the prefix to the front of every entry.
#
function(tribits_list_to_string LIST PREFIX OUTPUT_STRING)
  set(LIST_STRING "")

  foreach(ITEM ${LIST})
    set(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
  endforeach()

  set(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
endfunction()

#  This function will take a list of libraries and turn it into a space
#  separated string. In this case though the prefix is not always added
#  to the front of each entry as libraries can be specified either as a
#  name of a library to find or the absolute path to the library file
#  with any decorations the system uses. When an absolute path is given
#  the entry is used verbatim.
#
function(tribits_library_list_to_string LIST PREFIX OUTPUT_STRING)
  set(LIST_STRING "")

  foreach(ITEM ${LIST})
    string(SUBSTRING ${ITEM} 0 1 OPTION_FLAG)
    if(EXISTS ${ITEM} OR OPTION_FLAG STREQUAL "-")
      set(LIST_STRING "${LIST_STRING} ${ITEM}")
    else()
      set(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
    endif()
  endforeach()

  set(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_write_flexible_package_client_export_files()
#
# Utility function for writing ``${PACKAGE_NAME}Config.cmake`` files for
# package ``${PACKAGE_NAME}`` with some greater flexibility than what is
# provided by the function ``tribits_write_package_client_export_files()`` and
# to allow unit testing the generation of these files..
#
# Usage::
#
#   tribits_write_flexible_package_client_export_files(
#     PACKAGE_NAME <packageName>
#     [EXPORT_FILE_VAR_PREFIX <exportFileVarPrefix>]
#     [PACKAGE_CONFIG_FOR_BUILD_BASE_DIR <packageConfigForBuildBaseDir>]
#     [PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR <packageConfigForInstallBaseDir>]
#     )
#
# The arguments are:
#
#   ``PACKAGE_NAME <packageName>``
#
#     Gives the name of the TriBITS package for which the export files should
#     be created.
#
#   ``EXPORT_FILE_VAR_PREFIX <exportFileVarPrefix>``
#
#     If specified, then all of the variables in the generated export files
#     will be prefixed with ``<exportFileVarPrefix>_`` instead of
#     ``<packageName>_``.
#
#   ``PACKAGE_CONFIG_FOR_BUILD_BASE_DIR <packageConfigForBuildBaseDir>``
#
#     If specified, then the package's ``<packageName>Config.cmake`` file and
#     supporting files will be written under the directory
#     ``<packageConfigForBuildBaseDir>/`` (and any subdirs that does exist
#     will be created).  The generated file ``<packageName>Config.cmake`` is
#     for usage of the package in the build tree (not the install tree) and
#     points to include directories and libraries in the build tree.
#
#   ``PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR <packageConfigForInstallBaseDir>``
#
#     If specified, then the package's ``<packageName>Config_install.cmake``
#     file and supporting files will be written under the directory
#     ``<packageConfigForInstallBaseDir>/`` (and any subdirs that does exist
#     will be created).  The file ``${PACKAGE_NAME}Config_install.cmake`` is
#     meant to be installed renamed as ``<packageName>Config.cmake`` in the
#     install tree and it points to installed include directories and
#     libraries.
#
# NOTE: This function does *not* contain any ``install()`` command itself
# because CMake will not allow those to even be present in scripting mode that
# is used for unit testing this function.  Instead, the commands to install
# the files are added by the function
# ``tribits_write_package_client_export_files_install_targets()``.
#
function(tribits_write_flexible_package_client_export_files)

  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    message("\ntribits_write_flexible_package_client_export_files(${ARGN})")
  endif()

  #
  # A) Process the command-line arguments
  #

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     "WRITE_INSTALL_CMAKE_CONFIG_FILE"
     #one_value_keywords
     "PACKAGE_NAME;EXPORT_FILE_VAR_PREFIX;PACKAGE_CONFIG_FOR_BUILD_BASE_DIR;PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR"
     #multi_value_keywords
     ""
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  if (NOT ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES)
    message(SEND_ERROR "Error: Can't generate export dependency files because"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES is not ON!")
    return()
  endif()

  set(PACKAGE_NAME ${PARSE_PACKAGE_NAME})
  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    print_var(PACKAGE_NAME)
  endif()

  set(EXPORT_FILE_VAR_PREFIX ${PACKAGE_NAME})
  if (PARSE_EXPORT_FILE_VAR_PREFIX)
    set(EXPORT_FILE_VAR_PREFIX ${PARSE_EXPORT_FILE_VAR_PREFIX})
  endif()
  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    print_var(EXPORT_FILE_VAR_PREFIX)
  endif()

  #
  # B) Get the set of upstream packages for this package that are enabled,
  # libraries, library dirs, and include dirs
  #

  set(FULL_PACKAGE_SET "")
  set(FULL_LIBRARY_SET "")

  set(SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM TRUE)
  if (${PACKAGE_NAME}_INCLUDE_DIRS)
    set(FULL_INCLUDE_DIRS_SET ${${PACKAGE_NAME}_INCLUDE_DIRS})
    set(FULL_LIBRARY_DIRS_SET ${${PACKAGE_NAME}_LIBRARY_DIRS})
    set(SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM FALSE)
  else()
    set(FULL_INCLUDE_DIRS_SET "")
    set(FULL_LIBRARY_DIRS_SET "")
  endif()

  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    print_var(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES)
  endif()

  foreach(TRIBITS_PACKAGE ${${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES})

    if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
      print_var(TRIBITS_PACKAGE)
      if (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
        print_var(${TRIBITS_PACKAGE}_HAS_NATIVE_LIBRARIES_TO_INSTALL)
      endif()
    endif()

    set(APPEND_THE_PACKAGE TRUE)
    set(APPEND_THE_PACKAGE_LIBS TRUE)

    if (NOT ${TRIBITS_PACKAGE}_HAS_NATIVE_LIBRARIES_TO_INSTALL)
      set(APPEND_THE_PACKAGE_LIBS FALSE)
    endif()

    if (APPEND_THE_PACKAGE)
      list(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      if (APPEND_THE_PACKAGE_LIBS)
        append_set(FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
        if (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
          append_set(FULL_INCLUDE_DIRS_SET ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
          append_set(FULL_LIBRARY_DIRS_SET ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
        endif()
      else()
        if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
          message("-- " "Skipping adding the package libs!")
        endif()
      endif()
    else()
      if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
        message("-- " "Skipping adding the package!")
      endif()
    endif()

    if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
      print_var(FULL_PACKAGE_SET)
      print_var(FULL_LIBRARY_SET)
      if (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
        print_var(FULL_INCLUDE_DIRS_SET)
        print_var(FULL_LIBRARY_DIRS_SET)
      endif()
    endif()

  endforeach()

  # Must prepend the current package and its libraries itself so that we get
  # its TPLs libraries. However, if the current package has no native
  # libraries (yet), then there is no point in listing the package or its
  # TPLs.  Why would a package list TPLs (with actual libraries) if itself
  # does not have libraries to export?  Note, this does not affect internal
  # tests and examples which could have TPLs but no native libraries.
  if (${PACKAGE_NAME}_LIBRARIES AND ${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES_TO_INSTALL)
    prepend_set(FULL_PACKAGE_SET ${PACKAGE_NAME})
    prepend_set(FULL_LIBRARY_SET ${${PACKAGE_NAME}_LIBRARIES})
  endif()

  if (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
     if (FULL_INCLUDE_DIRS_SET)
       list(REMOVE_DUPLICATES FULL_INCLUDE_DIRS_SET)
     endif()
     if (FULL_LIBRARY_DIRS_SET)
       list(REMOVE_DUPLICATES FULL_LIBRARY_DIRS_SET)
     endif()
  endif()

  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    message("-- " "*** Final sets of packages, libs, include dirs, and lib dirs:")
    print_var(FULL_PACKAGE_SET)
    print_var(FULL_LIBRARY_SET)
    print_var(FULL_INCLUDE_DIRS_SET)
    print_var(FULL_LIBRARY_DIRS_SET)
  endif()

  #
  # C) Get the set of TPLs for this package that are enabled
  #

  # C.1) Get the set of enabled TPLs

  set(FULL_TPL_SET "")
  foreach(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    list(APPEND FULL_TPL_SET ${${TRIBITS_PACKAGE}_LIB_REQUIRED_DEP_TPLS})
    set(OPTIONAL_TPLS ${${TRIBITS_PACKAGE}_LIB_OPTIONAL_DEP_TPLS})
    foreach(TPL ${OPTIONAL_TPLS})
      # Only add if support for the optional TPL is enabled in this
      # package.  Don't just check if the TPL is enabled!
      if(${TRIBITS_PACKAGE}_ENABLE_${TPL})
        list(APPEND FULL_TPL_SET ${TPL})
      endif()
    endforeach()
  endforeach()
  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    print_var(FULL_TPL_SET)
  endif()

  # C.2) Sort the TPLs according to the master TPL list

  #We will use the complete list of supported tpls for the project
  #to help us create a properly ordered list of tpls.
  if (FULL_TPL_SET)
    set(ORDERED_FULL_TPL_SET ${FULL_TPL_SET})
    tribits_sort_list_according_to_master_list("${${PROJECT_NAME}_REVERSE_TPLS}"
      ORDERED_FULL_TPL_SET)
  endif()

  if (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    print_var(ORDERED_FULL_TPL_SET)
  endif()

  #
  # D) Get the libraries, library dirs, and the include dirs for the
  # upstream enabled TPLs
  #

  set(${PACKAGE_NAME}_TPL_LIBRARIES "")
  set(${PACKAGE_NAME}_TPL_INCLUDE_DIRS "")
  set(${PACKAGE_NAME}_TPL_LIBRARY_DIRS "")
  foreach(TPL ${ORDERED_FULL_TPL_SET})
    list(APPEND ${PACKAGE_NAME}_TPL_LIBRARIES ${TPL_${TPL}_LIBRARIES})
    list(APPEND ${PACKAGE_NAME}_TPL_INCLUDE_DIRS ${TPL_${TPL}_INCLUDE_DIRS})
    list(APPEND ${PACKAGE_NAME}_TPL_LIBRARY_DIRS ${TPL_${TPL}_LIBRARY_DIRS})
  endforeach()

  # Generate a note discouraging editing of the <package>Config.cmake file
  set(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  #
  # E) Deal with the library rpath issues with shared libs
  #

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.

  if(BUILD_SHARED_LIBS)
    string(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${FULL_LIBRARY_DIRS_SET}")
    set(SHARED_LIB_RPATH_COMMAND
      ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  endif()

  #
  # F) Create the contents of the <Package>Config.cmake file for the build tree
  #

  tribits_generate_package_config_file_for_build_tree(${PACKAGE_NAME})

  #
  # G) Create <Package>Config_install.cmake file for the install tree
  #

  tribits_generate_package_config_file_for_install_tree(${PACKAGE_NAME})

endfunction()


# @FUNCTION: tribits_generate_package_config_file_for_build_tree()
#
# Called from tribits_write_flexible_package_client_export_files() to finish
# up generating text for and writing the file `<Package>Config.cmake` for the
# build tree.
#
# Usage::
#
#   tribits_generate_package_config_file_for_build_tree(<packageName>)
#
# These files get placed under <buildDir>/cmake_packages/<packageName>/
#
# That makes them easy to find by find_package() by adding
# <buildDir>/cmake_packages/ to CMAKE_PREFIX_PATH.
#
function(tribits_generate_package_config_file_for_build_tree  packageName)

  set(buildDirExtPkgsDir
     "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}")
  set(buildDirCMakePkgsDir
     "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_CMAKE_PKGS_DIR}")

  if (PARSE_PACKAGE_CONFIG_FOR_BUILD_BASE_DIR
      OR PARSE_PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR
    )
    # Custom code in configuration file (gets pulled from by configure_file()
    # below)
    set(PACKAGE_CONFIG_CODE "")

    tribits_append_dependent_package_config_file_includes(${packageName}
      EXT_PKG_CONFIG_FILE_BASE_DIR "${buildDirExtPkgsDir}"
      PKG_CONFIG_FILE_BASE_DIR "${buildDirCMakePkgsDir}"
      CONFIG_FILE_STR_INOUT PACKAGE_CONFIG_CODE )

    # Import build tree targets into applications.
    #
    # BMA: Export only the immediate libraries of this project to the
    # build tree. Should manage more carefully, checking that they are
    # targets of this project and not other libs.  Also, should
    # consider more careful recursive management of targets when there
    # are sub-packages.  We'd like to export per-package, but deps
    # won't be satisfied, so we export one file for the project for
    # now...
    if (PARSE_PACKAGE_CONFIG_FOR_BUILD_BASE_DIR)
      tribits_get_package_config_build_dir_targets_file(${packageName}
        "${PACKAGE_CONFIG_FOR_BUILD_BASE_DIR}" packageConfigBuildDirTargetsFile )
      string(APPEND PACKAGE_CONFIG_CODE
        "\n# Import ${packageName} targets\n"
        "include(\"${packageConfigBuildDirTargetsFile}\")")
    endif()

    tribits_set_compiler_vars_for_config_file(BUILD_DIR)

    if ("${CMAKE_CXX_FLAGS}" STREQUAL "")
      set(CMAKE_CXX_FLAGS_ESCAPED "")
    else()
      # Replace " by \".
      string(REGEX REPLACE "\"" "\\\\\"" CMAKE_CXX_FLAGS_ESCAPED ${CMAKE_CXX_FLAGS})
    endif()
    set(tribitsConfigFilesDir
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}")
    configure_file(
      "${tribitsConfigFilesDir}/TribitsPackageConfigTemplate.cmake.in"
      "${PARSE_PACKAGE_CONFIG_FOR_BUILD_BASE_DIR}/${packageName}Config.cmake"
      )
  endif()

endfunction()


# @FUNCTION: tribits_generate_package_config_file_for_install_tree()
#
# Called from tribits_write_flexible_package_client_export_files() to finish
# up generating text for and writing the file `<Package>Config_install.cmake`
# that will get installed.
#
# Usage::
#
#   tribits_generate_package_config_file_for_install_tree(<packageName>)
#
# The export files are typically installed in
#     <install-dir>/<lib-path>/cmake/<package-name>/.
#
# This file isn't generally useful inside the build tree so it is being
# "hidden" in the CMakeFiles directory.  It gets installed in a separate
# function.  (NOTE: The install target is added in a different function to
# allow this function to be unit tested in a cmake -P script.)
#
function(tribits_generate_package_config_file_for_install_tree  packageName)

  # Set the include and library directories relative to the location
  # at which the ${PROJECT_NAME}Config.cmake file is going to be
  # installed. Note the variable reference below is escaped so it
  # won't be replaced until a client project attempts to locate
  # directories using the installed config file. This is to deal with
  # installers that allow relocation of the install tree at *install*
  # time.
  # The export files are typically installed in
  #     <install-dir>/<lib-path>/cmake/<package-name>/.
  # The relative path to the installation dir is hence k*(../) + ../../, where
  # k is the number of components in <lib-path>. Extract those here.
  # This doesn't work if ${${PROJECT_NAME}_INSTALL_LIB_DIR} contains "./" or
  # "../" components, but really, it never did. All of this should actually be
  # handled by CMake's configure_package_config_file().
  string(REPLACE "/" ";" PATH_LIST ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  set(RELATIVE_PATH "../..")
  foreach(PATH ${PATH_LIST})
    set(RELATIVE_PATH "${RELATIVE_PATH}/..")
  endforeach()
  set(FULL_LIBRARY_DIRS_SET
    "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")
  set(FULL_INCLUDE_DIRS_SET
    "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")

  # Custom code in configuration file.
  set(PACKAGE_CONFIG_CODE "")

  tribits_append_dependent_package_config_file_includes(${packageName}
    EXT_PKG_CONFIG_FILE_BASE_DIR
      "\${CMAKE_CURRENT_LIST_DIR}/../../${${PROJECT_NAME}_BUILD_DIR_EXTERNAL_PKGS_DIR}"
    PKG_CONFIG_FILE_BASE_DIR "\${CMAKE_CURRENT_LIST_DIR}/.."
    CONFIG_FILE_STR_INOUT PACKAGE_CONFIG_CODE )

  # Import install targets
  string(APPEND PACKAGE_CONFIG_CODE
    "\n# Import ${packageName} targets\n"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/${packageName}Targets.cmake\")")

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.
  if (BUILD_SHARED_LIBS)
    set(SHARED_LIB_RPATH_COMMAND
      ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR}
      )
  endif()

  tribits_set_compiler_vars_for_config_file(INSTALL_DIR)

  if (PARSE_PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR)
    configure_file(
      "${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in"
      "${PARSE_PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR}/${packageName}Config_install.cmake"
      )
  endif()

endfunction()


# @FUNCTION: tribits_append_dependent_package_config_file_includes()
#
# Append the includes for upstream external packages (TPLs) and internal
# packages to a `<Package>Config.cmake` file string.
#
# Usage::
#
#   tribits_append_dependent_package_config_file_includes(
#     <packageName>
#     EXT_PKG_CONFIG_FILE_BASE_DIR <extPkgconfigFileBaseDir>
#     PKG_CONFIG_FILE_BASE_DIR <pkgConfigFileBaseDir>
#     CONFIG_FILE_STR_INOUT <configFileStrInOut>
#     )
#
function(tribits_append_dependent_package_config_file_includes packageName)

  # Parse input

  cmake_parse_arguments(
     PARSE  #prefix
     ""  #options
     #one_value_keywords
     "EXT_PKG_CONFIG_FILE_BASE_DIR;PKG_CONFIG_FILE_BASE_DIR;CONFIG_FILE_STR_INOUT"
     "" #multi_value_keywords
     ${ARGN}
     )
  tribits_check_for_unparsed_arguments()

  set(extPkgConfigFileBaseDir "${PARSE_EXT_PKG_CONFIG_FILE_BASE_DIR}")
  set(pkgConfigFileBaseDir "${PARSE_PKG_CONFIG_FILE_BASE_DIR}")
  set(configFileStr "${${PARSE_CONFIG_FILE_STR_INOUT}}")

  # Include configurations of dependent packages
  string(APPEND configFileStr
    "# Include configuration of dependent packages\n")
  foreach(depPkg IN LISTS ${packageName}_FULL_ENABLED_DEP_PACKAGES)
    set(cmakePkgDir "${pkgConfigFileBaseDir}/${depPkg}")
    string(APPEND configFileStr
      "include(\"${cmakePkgDir}/${depPkg}Config.cmake\")\n")
  endforeach()

  # Include configurations of dependent external packages/TPLs
  string(APPEND configFileStr
    "\n# Include configuration of dependent external packages/TPls\n")
  foreach(depTpl IN LISTS ${packageName}_LIB_REQUIRED_DEP_TPLS)
    if (TARGET ${depTpl}::all_libs)
      set(cmakeTplDir "${extPkgConfigFileBaseDir}/${depTpl}")
      string(APPEND configFileStr
        "include(\"${cmakeTplDir}/${depTpl}Config.cmake\")\n")
    endif()
  endforeach()
  foreach(depTpl IN LISTS ${packageName}_LIB_OPTIONAL_DEP_TPLS)
    if (${packageName}_ENABLE_${depTpl} AND TARGET ${depTpl}::all_libs)
      set(cmakeTplDir "${extPkgConfigFileBaseDir}/${depTpl}")
      string(APPEND configFileStr
        "include(\"${cmakeTplDir}/${depTpl}Config.cmake\")\n")
    endif()
  endforeach()
  # NOTE: Above, every TPL does not have a <tplName>Config.cmake file written
  # for it.  For example, special TPLs like "MPI" don't have this file created
  # or have an MPI::all_libs target corrected.  Therefore, we check for the
  # defintion <tplName>::all_libs before we include the file above.

  # Set the output
  set(${PARSE_CONFIG_FILE_STR_INOUT} "${configFileStr}" PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_write_package_client_export_files_install_targets()
#
# Create the ``<Package>ConfigTargets.cmake`` file and install rules and the
# install() target for the previously generated
# ``<Package>Config_install.cmake`` files generated by the
# `tribits_write_flexible_package_client_export_files()`_ function.
#
# Usage::
#
#   tribits_write_package_client_export_files_install_targets(
#     PACKAGE_NAME <packageName>
#     PACKAGE_CONFIG_FOR_BUILD_BASE_DIR <packageConfigForBuildBaseDir>
#     PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR <packageConfigForInstallBaseDir>
#     )
#
# The install() commands must be in a different subroutine or CMake will not
# allow you to call the routine, even if you if() it out!
#
function(tribits_write_package_client_export_files_install_targets)

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     ""
     #one_value_keywords
     "PACKAGE_NAME;PACKAGE_CONFIG_FOR_BUILD_BASE_DIR;PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR"
     #multi_value_keywords
     ""
     ${ARGN}
     )

  set(PACKAGE_NAME ${PARSE_PACKAGE_NAME})

  if (PARSE_PACKAGE_CONFIG_FOR_BUILD_BASE_DIR)
    tribits_get_package_config_build_dir_targets_file(${PACKAGE_NAME}
      "${PARSE_PACKAGE_CONFIG_FOR_BUILD_BASE_DIR}" packageConfigBuildDirTargetsFile )
    export(
      EXPORT ${PACKAGE_NAME}
      NAMESPACE ${PACKAGE_NAME}::
      FILE "${packageConfigBuildDirTargetsFile}" )
  endif()

  if (PARSE_PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR)
    install(
      FILES
        "${PARSE_PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR}/${PACKAGE_NAME}Config_install.cmake"
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
      RENAME ${PACKAGE_NAME}Config.cmake
      )
    install(
      EXPORT ${PACKAGE_NAME}
      NAMESPACE ${PACKAGE_NAME}::
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
      FILE "${PACKAGE_NAME}Targets.cmake" )
  endif()

endfunction()


# Function to return the full path the targets file for the
# <Package>Config.cmake file in the build tree.
#
function(tribits_get_package_config_build_dir_targets_file  PACKAGE_NAME
    packageConfigForBuildBaseDir  packageConfigBuildDirTargetsFileOut
  )
  set(${packageConfigBuildDirTargetsFileOut}
    "${PACKAGE_CONFIG_FOR_BUILD_BASE_DIR}/${PACKAGE_NAME}Targets.cmake"
    PARENT_SCOPE )
endfunction()


# Generate the ${PACKAGE_NAME}Config.cmake file for package PACKAGE_NAME.
#
# ToDo: Finish documentation!
#
function(tribits_write_package_client_export_files PACKAGE_NAME)

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES: ${PACKAGE_NAME}")
  endif()

  set(buildDirCMakePkgsDir
     "${${PROJECT_NAME}_BINARY_DIR}/${${PROJECT_NAME}_BUILD_DIR_CMAKE_PKGS_DIR}")

  set(EXPORT_FILES_ARGS PACKAGE_NAME ${PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("For package ${PACKAGE_NAME} creating ${PACKAGE_NAME}Config.cmake")
    endif()
    set(PACKAGE_CONFIG_FOR_BUILD_BASE_DIR
      "${buildDirCMakePkgsDir}/${PACKAGE_NAME}" )
    set(PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR
      "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles" )
    append_set(EXPORT_FILES_ARGS
      PACKAGE_CONFIG_FOR_BUILD_BASE_DIR "${PACKAGE_CONFIG_FOR_BUILD_BASE_DIR}"
      PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR "${PACKAGE_CONFIG_FOR_INSTALL_BASE_DIR}"
      )
  endif()

  tribits_write_flexible_package_client_export_files(${EXPORT_FILES_ARGS})

  tribits_write_package_client_export_files_install_targets(${EXPORT_FILES_ARGS})

endfunction()


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
  set(PACKAGE_LIST ${${PROJECT_NAME}_SE_PACKAGES})
  if (PACKAGE_LIST)
    list(REVERSE PACKAGE_LIST)
  endif()

  # Loop over all packages to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  set(FULL_PACKAGE_SET "")
  set(FULL_LIBRARY_SET "")
  set(FULL_INCLUDE_DIRS_SET "")
  set(FULL_LIBRARY_DIRS_SET "")
  foreach(TRIBITS_PACKAGE ${PACKAGE_LIST})
    if(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      list(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      list(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
      list(APPEND FULL_INCLUDE_DIRS_SET ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
      list(APPEND FULL_LIBRARY_DIRS_SET ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
    endif()
  endforeach()

  set(${PROJECT_NAME}_CONFIG_LIBRARIES ${FULL_LIBRARY_SET})

  # Reversing the tpl list so that the list of tpls will be produced in
  # order of most dependent to least dependent.
  if (${PROJECT_NAME}_TPLS)
    set(TPL_LIST ${${PROJECT_NAME}_TPLS})
    list(REVERSE TPL_LIST)
  endif()

  # Loop over all TPLs to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  set(FULL_TPL_SET "")
  set(FULL_TPL_LIBRARY_SET "")
  set(FULL_TPL_INCLUDE_DIRS_SET "")
  set(FULL_TPL_LIBRARY_DIRS_SET "")
  foreach(TPL ${TPL_LIST})
    if(TPL_ENABLE_${TPL})
      list(APPEND FULL_TPL_SET ${TPL})
      list(APPEND FULL_TPL_LIBRARY_SET ${TPL_${TPL}_LIBRARIES})
      list(APPEND FULL_TPL_INCLUDE_DIRS_SET ${TPL_${TPL}_INCLUDE_DIRS})
      list(APPEND FULL_TPL_LIBRARY_DIRS_SET ${TPL_${TPL}_LIBRARY_DIRS})
    endif()
  endforeach()

  # it is possible that tpls are in the same directory, to keep from
  # having a very long include path or library path we will strip out
  # any duplicates. This shouldn't affect which include or library is
  # found since the first instance of any path will be the one that is
  # kept.
  list(REMOVE_DUPLICATES FULL_TPL_INCLUDE_DIRS_SET)
  list(REMOVE_DUPLICATES FULL_TPL_LIBRARY_DIRS_SET)

  set(${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS ${FULL_TPL_INCLUDE_DIRS_SET})
  set(${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS ${FULL_TPL_LIBRARY_DIRS_SET})
  set(${PROJECT_NAME}_CONFIG_TPL_LIBRARIES ${FULL_TPL_LIBRARY_SET})

  #
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree
  # and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  set(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  # Config file for setting variables and finding include/library paths from
  # the build directory
  set(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  set(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    string(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND
     "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
    set(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  endif()

  # Custom code in configuration file.
  set(PROJECT_CONFIG_CODE "")

  #  # Export targets from the build tree.
  #  if(FULL_LIBRARY_SET)
  #    list(SORT FULL_LIBRARY_SET)
  #    list(REMOVE_DUPLICATES FULL_LIBRARY_SET)
  #    set(FULL_LIBRARY_TARGET_SET)
  #    foreach(LIB_ELE ${FULL_LIBRARY_SET})
  #      if (TARGET ${LIB_ELE})
  #        list(APPEND FULL_LIBRARY_TARGET_SET ${LIB_ELE})
  #      endif()
  #    endforeach()
  #    export(TARGETS ${FULL_LIBRARY_TARGET_SET} FILE
  #      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
  #    # Import the targets in applications.
  #    set(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
  ## Import ${PROJECT_NAME} targets
  #if(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  #  set(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  #  include(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
  #endif()
  #")
  #  endif()

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
  set(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS
    "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
  set(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS
    "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")

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



macro(tribits_set_compiler_var_for_config_file LANG FOR_DIR)
  if (NOT "${CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE_${FOR_DIR}}" STREQUAL "")
    set(CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE
      "${CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE_${FOR_DIR}}")
  else()
    set(CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE
      "${CMAKE_${LANG}_COMPILER}")
  endif()
  #message("${FOR_DIR}: CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE='${CMAKE_${LANG}_COMPILER_FOR_CONFIG_FILE}'")
endmacro()


macro(tribits_set_compiler_vars_for_config_file FOR_DIR)
  tribits_set_compiler_var_for_config_file(CXX ${FOR_DIR})
  tribits_set_compiler_var_for_config_file(C ${FOR_DIR})
  tribits_set_compiler_var_for_config_file(Fortran ${FOR_DIR})
endmacro()
