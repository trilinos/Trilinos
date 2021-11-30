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

#
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

#
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

#
# CMAKE_CURRENT_LIST_DIR is not defined in CMake versions < 2.8.3, but the
# Trilinos writes paths that use the value of that variable to this file.
# Make sure it is available at *find_package* time. Note that all variable
# references in the code snippet are escaped. This is to keep them from
# being evaluated until they are actually in the install tree. This is
# done to handle movable install trees.
#
# This function defines the variable
# DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_CODE_SNIPPET in the caller's scope
# as a string that can be referenced from CONFIGURE_FILE input files
# to ensure that the CMAKE_CURRENT_LIST_DIR will be defined on the installation
# target machine, even if it has an older version of cmake.
#
function(tribits_set_define_cmake_current_list_dir_code_snippet)
  set(DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET "
# Include guard
if (${EXPORT_FILE_VAR_PREFIX}_CONFIG_INCLUDED)
  return()
endif()
set(${EXPORT_FILE_VAR_PREFIX}_CONFIG_INCLUDED TRUE)

# Make sure CMAKE_CURRENT_LIST_DIR is usable
if (NOT DEFINED CMAKE_CURRENT_LIST_DIR)
  get_filename_component(_THIS_SCRIPT_PATH \${CMAKE_CURRENT_LIST_FILE} PATH)
  set(CMAKE_CURRENT_LIST_DIR \${_THIS_SCRIPT_PATH})
endif()
"
  PARENT_SCOPE )
endfunction()


#
# @FUNCTION: tribits_write_flexible_package_client_export_files()
#
# Utility function for writing ``${PACKAGE_NAME}Config.cmake`` and/or the
# ``Makefile.export.${PACKAGE_NAME}`` files for package ``${PACKAGE_NAME}``
# with some greater flexibility than what is provided by the function
# ``tribits_write_package_client_export_files()``.
#
# Usage::
#
#   tribits_write_flexible_package_client_export_files(
#     PACKAGE_NAME <packageName>
#     [EXPORT_FILE_VAR_PREFIX <exportFileVarPrefix>]
#     [WRITE_CMAKE_CONFIG_FILE <cmakeConfigFileFullPath>]
#     [WRITE_EXPORT_MAKEFILE <exportMakefileFileFullPath>]
#     [WRITE_INSTALL_CMAKE_CONFIG_FILE]
#     [WRITE_INSTALL_EXPORT_MAKEFILE]
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
#   ``WRITE_CMAKE_CONFIG_FILE <cmakeConfigFileFullPath>``
#
#     If specified, then the package's (``<packageName>``) cmake configure
#     export file for use by external CMake client projects will be created as
#     the file ``<cmakeConfigFileFullPath>``.  NOTE: the argument should be
#     the full path!
#
#   ``WRITE_EXPORT_MAKEFILE <exportMakefileFileFullPath>``
#
#     If specified, then the package's (``<packageName>``) export makefile for
#     use by external Makefile client projects will be created in the file
#     <exportMakefileFileFullPath>.  NOTE: the argument should be the full
#     path!
#
#   ``WRITE_INSTALL_CMAKE_CONFIG_FILE``
#
#     If specified, then the package's (``<packageName>``) install cmake
#     configured export file will be installed in to the install tree as well.
#     The name and location of this file is hard-coded.
#
#   ``WRITE_INSTALL_EXPORT_MAKEFILE``
#
#     If specified, then the package's (``<packageName>``) install export
#     makefile to be installed into the install tree as well.  The name and
#     location of this file is hard-coded.
#
# NOTE: The arguments to this function may look strange but the motivation is
# to support very specialized use cases such as when a TriBITS package needs
# to generate an export makefile for a given package but the name of the
# export makefile must be different and use different variable name prefixes.
# The particular use case is when wrapping an external autotools project that
# depends on Trilinos and needs to read in the ``Makefile.export.Trilinos``
# file but this file needs to be generated for a subset of enabled packages on
# the fly during a one-pass configure.
#
# NOTE: This function does *not* contain the ``install()`` commands because
# CMake will not allow those to even be present in scripting mode that is used
# for unit testing this function.  Instead, the files to be installed are only
# generated in the build tree and the install targets are added else where.
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
     "WRITE_INSTALL_CMAKE_CONFIG_FILE;WRITE_INSTALL_EXPORT_MAKEFILE"
     #one_value_keywords
     ""
     #multi_value_keywords
     "PACKAGE_NAME;WRITE_CMAKE_CONFIG_FILE;WRITE_EXPORT_MAKEFILE;EXPORT_FILE_VAR_PREFIX"
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

  tribits_set_define_cmake_current_list_dir_code_snippet()

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
  # F) Create the contents of the cmake Config.cmake file for the build tree
  #
  # Creating this file in the base dir of the package since it is possible
  # that the cmake returned path where the file was found would be useful for
  # a package, and having to dig through the hiding that is done wouldn't be
  # nice.
  #

  if (PARSE_WRITE_CMAKE_CONFIG_FILE)
    # Custom code in configuration file.
    set(PACKAGE_CONFIG_CODE "")

    # Include configurations of dependent packages
    foreach(DEP_PACKAGE ${${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES})
      # Could use file(RELATIVE_PATH ...), but probably not necessary
      # since unlike install trees, build trees need not be relocatable
      set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
include(\"${${DEP_PACKAGE}_BINARY_DIR}/${DEP_PACKAGE}Config.cmake\")"
        )
    endforeach()

    # Import build tree targets into applications.
    #
    # BMA: Export only the immediate libraries of this project to the
    # build tree. Should manage more carefully, checking that they are
    # targets of this project and not other libs.  Also, should
    # consider more careful recursive management of targets when there
    # are sub-packages.  We'd like to export per-package, but deps
    # won't be satisfied, so we export one file for the project for
    # now...
    if(${TRIBITS_PACKAGE}_HAS_NATIVE_LIBRARIES_TO_INSTALL)
      export(TARGETS ${${PACKAGE_NAME}_LIBRARIES} FILE
	"${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Targets.cmake" APPEND)
      set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PACKAGE_NAME} targets
include(\"${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")"
      )
    endif()

    tribits_set_compiler_vars_for_config_file(BUILD_DIR)

    if ("${CMAKE_CXX_FLAGS}" STREQUAL "")
      set(CMAKE_CXX_FLAGS_ESCAPED "")
    else()
      # Replace " by \".
      string(REGEX REPLACE "\"" "\\\\\"" CMAKE_CXX_FLAGS_ESCAPED ${CMAKE_CXX_FLAGS})
    endif()
    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in
      "${PARSE_WRITE_CMAKE_CONFIG_FILE}"
      )
  endif()

  #
  # G) Create the export makefile for the build tree
  #
  # This is the equivalent of the cmake version only slightly changed so that
  # it can be directly imported into a Makefile.
  #

  if(PARSE_WRITE_EXPORT_MAKEFILE)

    tribits_list_to_string("${FULL_LIBRARY_SET}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_FULL_LIBRARY_SET)
    tribits_list_to_string("${FULL_LIBRARY_DIRS_SET}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    tribits_list_to_string("${FULL_INCLUDE_DIRS_SET}" "-I" MAKEFILE_INCLUDE_DIRS)
    tribits_list_to_string("${${PACKAGE_NAME}_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PACKAGE_NAME}_TPL_INCLUDE_DIRS)
    tribits_list_to_string("${${PACKAGE_NAME}_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    tribits_library_list_to_string("${${PACKAGE_NAME}_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARIES)

    tribits_library_list_to_string("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    tribits_list_to_string("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    tribits_list_to_string("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")

    tribits_list_to_string("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    tribits_list_to_string("${ORDERED_FULL_TPL_SET}" "" MAKEFILE_ORDERED_FULL_TPL_SET)

    # create an upper case name of the package so that we can make deprecated
    # versions of them to help people transistioning from the autotools
    # version diagnose any missed variables.
    string(TOUPPER ${EXPORT_FILE_VAR_PREFIX} EXPORT_FILE_VAR_PREFIX_UPPER)

    assert_defined(${PROJECT_NAME}_TRIBITS_DIR)
    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in
      "${PARSE_WRITE_EXPORT_MAKEFILE}"
      )
  endif()

  #
  # H) Create the cmake Config.cmake file for the install tree.
  #

  # This file isn't generally useful inside the build tree so it is being
  # "hidden" in the CMakeFiles directory. It will be placed in the base
  # install directory for ${PROJECT_NAME} when installed.

  # NOTE: The install target is added in a different function to allow this
  # function to be unit tested in a cmake -P script.

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
  set(FULL_LIBRARY_DIRS_SET "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")
  set(FULL_INCLUDE_DIRS_SET "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")

  # Custom code in configuration file.
  set(PACKAGE_CONFIG_CODE "")

  if (${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES)
    set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Include configuration of dependent packages")
  endif()
  foreach(DEP_PACKAGE ${${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES})
    set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
include(\"\${CMAKE_CURRENT_LIST_DIR}/../${DEP_PACKAGE}/${DEP_PACKAGE}Config.cmake\")"
)
  endforeach()
  if(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES)
    set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}\n")
  endif()

  # Import install tree targets into applications.
  get_property(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS)
  if(HAS_INSTALL_TARGETS)
    set(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PACKAGE_NAME} targets
include(\"\${CMAKE_CURRENT_LIST_DIR}/${PACKAGE_NAME}Targets.cmake\")"
)
  endif()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    set(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()

  tribits_set_compiler_vars_for_config_file(INSTALL_DIR)

  if (PARSE_WRITE_INSTALL_CMAKE_CONFIG_FILE)

    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
      )

  endif()

  #
  # I) Write the export makefile for the install tree
  #

  if (PARSE_WRITE_INSTALL_EXPORT_MAKEFILE)

    # Generated Make imports must use CMAKE_INSTALL_PREFIX, rather
    # than the more platform friendly method of locating the libraries
    # and includes using the config file path above. The underlying
    # assumption here is that a generator that uses
    # CMAKE_INSTALL_PREFIX is being used.
    set(FULL_LIBRARY_DIRS_SET ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
    set(FULL_INCLUDE_DIRS_SET ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

    tribits_list_to_string("${FULL_LIBRARY_DIRS_SET}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    tribits_list_to_string("${FULL_INCLUDE_DIRS_SET}" "-I" MAKEFILE_INCLUDE_DIRS)

    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      )
  endif()

endfunction()


#
# Set the install targets for the package config and export makefiles.
#
# The install() commands must be in a different subroutine or CMake will not
# allow you to call the routine, even if you if() it out!
#

function(tribits_write_project_client_export_files_install_targets PACKAGE_NAME)

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
      RENAME ${PACKAGE_NAME}Config.cmake
      )

    if(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES_TO_INSTALL)
      install(
        EXPORT ${PACKAGE_NAME}
        DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
        FILE ${PACKAGE_NAME}Targets.cmake
        )
    endif()
  endif()

  if(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)

    install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PACKAGE_NAME}
      )
  endif()

endfunction()


#
# Generate the ${PACKAGE_NAME}Config.cmake and/or the Makefile.export.${PACKAGE_NAME}
# for package PACKAGE_NAME.
#
# ToDo: Finish documentation!
#

function(tribits_write_package_client_export_files PACKAGE_NAME)

  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nTRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES: ${PACKAGE_NAME}")
  endif()

  set(EXPORT_FILES_ARGS PACKAGE_NAME ${PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    set(WRITE_CMAKE_CONFIG_FILE
      ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake)
    if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("For package ${PACKAGE_NAME} creating ${WRITE_CMAKE_CONFIG_FILE}")
    endif()
    append_set(EXPORT_FILES_ARGS
      WRITE_CMAKE_CONFIG_FILE "${WRITE_CMAKE_CONFIG_FILE}"
      WRITE_INSTALL_CMAKE_CONFIG_FILE)
  endif()

  if(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    set(WRITE_EXPORT_MAKEFILE
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PACKAGE_NAME})
    if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("For package ${PACKAGE_NAME} creating ${WRITE_EXPORT_MAKEFILE}")
    endif()
    append_set(EXPORT_FILES_ARGS
      WRITE_EXPORT_MAKEFILE "${WRITE_EXPORT_MAKEFILE}"
      WRITE_INSTALL_EXPORT_MAKEFILE)
  endif()

  tribits_write_flexible_package_client_export_files(${EXPORT_FILES_ARGS})

  tribits_write_project_client_export_files_install_targets(${PACKAGE_NAME})

endfunction()


#
# Write the outer TriBITS project configure and/or export makefiles
#
# If ${PROJECT_NAME}_VERSION is not set or is '' on input, then it will be set
# to 0.0.0 in order to create the ${PROJECT_NAME}ConfigVersion.cmake file.
#
# ToDo: Finish documentation!
#

function(tribits_write_project_client_export_files)

  set(EXPORT_FILE_VAR_PREFIX ${PROJECT_NAME})
  tribits_set_define_cmake_current_list_dir_code_snippet()

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
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  set(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  #Config file for setting variables and finding include/library paths from the build directory
  set(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  set(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    string(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
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
    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
      ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake )
  endif()

  if(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARIES)
    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)
    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS)
    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    tribits_library_list_to_string("${${PROJECT_NAME}_CONFIG_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARIES)

    tribits_library_list_to_string("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    tribits_list_to_string("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    tribits_list_to_string("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")

    tribits_list_to_string("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    tribits_list_to_string("${FULL_TPL_SET}" "" MAKEFILE_FULL_TPL_SET)

    if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
      # In TribitsProjectConfigTemplate.cmake.in, we would like to preserve
      # ${}-variables after the conversion to
      # TribitsProjectConfigTemplate.cmake. To this end, one typically uses the
      # @-syntax for variables. That doesn't support nested variables, however.
      # Use ${PDOLLAR} as a workaround, cf.
      # <http://www.cmake.org/pipermail/cmake/2013-April/054341.html>.
      set(PDOLLAR "$")
      configure_file(
        ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
        ${PROJECT_BINARY_DIR}/Makefile.export.${PROJECT_NAME})
    endif()
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
  set(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
  set(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS "\${CMAKE_CURRENT_LIST_DIR}/${RELATIVE_PATH}/${${PROJECT_NAME}_INSTALL_LIB_DIR}")

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  if(BUILD_SHARED_LIBS)
    set(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()

  # Custom code in configuration file.
  set(PROJECT_CONFIG_CODE "")

  tribits_set_compiler_vars_for_config_file(INSTALL_DIR)

  if (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
      ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

    install(
      FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      RENAME ${PROJECT_NAME}Config.cmake
      )
  endif()

  if(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the install tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    # Generated Make imports must use CMAKE_INSTALL_PREFIX, rather
    # than the more platform friendly method of locating the libraries
    # and includes using the config file path above. The underlying
    # assumption here is that a generator that uses
    # CMAKE_INSTALL_PREFIX is being used.
    set(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
    set(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})

    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    tribits_list_to_string("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)

    configure_file(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
      ${PROJECT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install )

    install(
      FILES ${PROJECT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PROJECT_NAME}
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
    FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
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
