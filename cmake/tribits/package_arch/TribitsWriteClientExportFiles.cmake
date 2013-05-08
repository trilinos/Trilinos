# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

INCLUDE(TribitsGeneralMacros)

###
### WARNING: See "NOTES TO DEVELOPERS" at the bottom of the file
### TribitsPackageMacros.cmake!
###

#
#  This function will take a list and turn it into a space separated string
#  adding the prefix to the front of every entry.
#
FUNCTION(TRIBITS_LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")

  FOREACH(ITEM ${LIST})
    SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

#
#  This function will take a list of libraries and turn it into a space
#  separated string. In this case though the prefix is not always added
#  to the front of each entry as libraries can be specified either as a
#  name of a library to find or the absolute path to the library file
#  with any decorations the system uses. When an absolute path is given
#  the entry is used verbatim.
#
FUNCTION(TRIBITS_LIBRARY_LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")

  FOREACH(ITEM ${LIST})
    STRING(SUBSTRING ${ITEM} 0 1 OPTION_FLAG)
    IF(EXISTS ${ITEM} OR OPTION_FLAG STREQUAL "-")
      SET(LIST_STRING "${LIST_STRING} ${ITEM}")
    ELSE()
      SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
    ENDIF()
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

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
FUNCTION(TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET)
  SET(DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET "
IF (NOT DEFINED CMAKE_CURRENT_LIST_DIR)
  GET_FILENAME_COMPONENT(_THIS_SCRIPT_PATH \${CMAKE_CURRENT_LIST_FILE} PATH)
  SET(CMAKE_CURRENT_LIST_DIR \${_THIS_SCRIPT_PATH})
ENDIF()
"
  PARENT_SCOPE )
ENDFUNCTION()


#
# @BEGIN(FUNCTION) TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES
#
# TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(
#   PACKAGE_NAME <pakageName>
#   [EXPORT_FILE_VAR_PREFIX <exportFileVarPrefix>]
#   [WRITE_CMAKE_CONFIG_FILE <cmakeConfigFileFullPath>]
#   [WRITE_EXPORT_MAKLEFILE <exportMakefileFileFullPath>]
#   [WRITE_INSTALL_CMAKE_CONFIG_FILE]
#   [WRITE_INSTALL_EXPORT_MAKLEFILE]
#   )
#
# Utility function for writing ${PACKAGE_NAME}Config.cmake and/or the
# Makefile.export.${PACKAGE_NAME} for package PACKAGE_NAME with some greater
# flexibility than TRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES()
#
# The arguments are:
#
# PACKAGE_NAME <pakageName>
#
#   Gives the name of the TriBITS package for which the export files should be
#   created for.
#
# EXPORT_FILE_VAR_PREFIX <exportFileVarPrefix>
#
#   If specified, then all of the varibles in the generated export files will
#   be prefixed with "<exportFileVarPrefix>_" instead of "${PACKAGE_NAME}_".
#   This is to provide flexibility.
#
# WRITE_CMAKE_CONFIG_FILE <cmakeConfigFileFullPath>
#
#   If specified, then the package <packageName>'s cmake configure export file
#   for extenral CMake client projects will be created in the file
#   <cmakeConfigFileFullPath>.  NOTE: the argument should be the full path!
#
# WRITE_EXPORT_MAKLEFILE <exportMakefileFileFullPath>
#
#   If specified, then the package <packageName>'s cmake configure export file
#   for external Makefile client projects will be created in the file
#   <exportMakefileFileFullPath>.  NOTE: the argument should be the full path!
#
# WRITE_INSTALL_CMAKE_CONFIG_FILE
#
#   If specified, then the package <packageName>'s install cmake configure
#   export to be installed will be written.  The name and location of this
#   file is hard-coded.
#
# WRITE_INSTALL_EXPORT_MAKLEFILE
#
#   If specified, then the package <packageName>'s install export makefile to
#   be installed will be written.  The name and location of this file is
#   hard-coded.
#
# NOTE: The arguments to this function may look strange but the motivation is
# to support versy speicalized use cases such as when a TriBITS package needs
# to generate an export makefile for a given package but name the export
# makefile differently and use different variable name prefixes.  The
# particular driver use case is when wrapping an external autotools project
# that depends on Trilinos and needs to read in the Makefile.export.Trilinos
# file but this file needs to be generated for a subset of enabled packages on
# the fly during a one-pass configure.
#
# NOTE: This function does *not* contain the the INSTALL() commands because
# CMake will not allow those to even be present in scripting mode that is used
# for unit testing this function.
#
# @END(FUNCTION) TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES
#

FUNCTION(TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES)

  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    MESSAGE("\nTRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(${ARGN})")
  ENDIF()

  TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET()

  #
  # A) Process the command-line arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "PACKAGE_NAME;WRITE_CMAKE_CONFIG_FILE;WRITE_EXPORT_MAKLEFILE;EXPORT_FILE_VAR_PREFIX"
     #options
     "WRITE_INSTALL_CMAKE_CONFIG_FILE;WRITE_INSTALL_EXPORT_MAKLEFILE"
     ${ARGN}
     )

  IF (NOT ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES)
    MESSAGE(SEND_ERROR "Error: Can't generate export depenency files because"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES is not ON!")
    RETURN()
  ENDIF()

  SET(PACKAGE_NAME ${PARSE_PACKAGE_NAME})
  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    PRINT_VAR(PACKAGE_NAME)
  ENDIF()

  SET(EXPORT_FILE_VAR_PREFIX ${PACKAGE_NAME})
  IF (PARSE_EXPORT_FILE_VAR_PREFIX)
    SET(EXPORT_FILE_VAR_PREFIX ${PARSE_EXPORT_FILE_VAR_PREFIX})
  ENDIF()
  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    PRINT_VAR(EXPORT_FILE_VAR_PREFIX)
  ENDIF()

  #
  # B) Get the set of upstream packages for this package that are enabled,
  # libraries, library dirs, and include dirs
  #

  SET(FULL_PACKAGE_SET "")
  SET(FULL_LIBRARY_SET "")

  SET(SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM TRUE)
  IF (${PACKAGE_NAME}_INCLUDE_DIRS)
    SET(FULL_INCLUDE_DIRS_SET ${${PACKAGE_NAME}_INCLUDE_DIRS})
    SET(FULL_LIBRARY_DIRS_SET ${${PACKAGE_NAME}_LIBRARY_DIRS})
    SET(SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM FALSE)
  ELSE()
    SET(FULL_INCLUDE_DIRS_SET "")
    SET(FULL_LIBRARY_DIRS_SET "")
  ENDIF()

  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    PRINT_VAR(${PACKAGE_NAME}_FULL_EXPORT_DEP_PACKAGES)
  ENDIF()

  FOREACH(TRIBITS_PACKAGE ${${PACKAGE_NAME}_FULL_EXPORT_DEP_PACKAGES})

    LIST(FIND ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES ${TRIBITS_PACKAGE}
      PACKAGE_IS_OPTIONAL_DEP)

    IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
      PRINT_VAR(TRIBITS_PACKAGE)
      PRINT_VAR(PACKAGE_IS_OPTIONAL_DEP)
      IF (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
        PRINT_VAR(${TRIBITS_PACKAGE}_HAS_NATIVE_LIBRARIES)
      ENDIF()
    ENDIF()

    SET(APPEND_THE_PACKAGE TRUE)
    SET(APPEND_THE_PACKAGE_LIBS TRUE)

    IF(PACKAGE_IS_OPTIONAL_DEP GREATER -1)
      IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
        PRINT_VAR(${PACKAGE_NAME}_ENABLE_${TRIBITS_PACKAGE})
      ENDIF()
      IF (NOT ${PACKAGE_NAME}_ENABLE_${TRIBITS_PACKAGE})
        SET(APPEND_THE_PACKAGE FALSE)
        SET(APPEND_THE_PACKAGE_LIBS FALSE)
      ENDIF()
    ENDIF()

    IF (NOT ${TRIBITS_PACKAGE}_HAS_NATIVE_LIBRARIES)
      SET(APPEND_THE_PACKAGE_LIBS FALSE)
    ENDIF()

    IF (APPEND_THE_PACKAGE)
      LIST(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      IF (APPEND_THE_PACKAGE_LIBS)
        APPEND_SET(FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
        IF (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
          APPEND_SET(FULL_INCLUDE_DIRS_SET ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
          APPEND_SET(FULL_LIBRARY_DIRS_SET ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
        ENDIF()
      ELSE()
        IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
          MESSAGE("-- " "Skipping adding the package libs!")
        ENDIF()
      ENDIF()
    ELSE()
      IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
        MESSAGE("-- " "Skipping adding the package!")
      ENDIF()
    ENDIF()

    IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
      PRINT_VAR(FULL_PACKAGE_SET)
      PRINT_VAR(FULL_LIBRARY_SET)
      IF (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
        PRINT_VAR(FULL_INCLUDE_DIRS_SET)
        PRINT_VAR(FULL_LIBRARY_DIRS_SET)
      ENDIF()
    ENDIF()

  ENDFOREACH()

  # Must prepend the current package and its libraries itself so that we get
  # its TPLs libraries. However, if the current pacakge has no native
  # libraries (yet), then there is no point in listing the package or its
  # TPLs.  Why would a package list TPLs (with actual libraries) if itself
  # does not have libraries to export?  Note, this does not affect internal
  # tests and examples which could have TPLs but no native libraries.
  IF (${PACKAGE_NAME}_LIBRARIES AND ${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES)
    PREPEND_SET(FULL_PACKAGE_SET ${PACKAGE_NAME})
    PREPEND_SET(FULL_LIBRARY_SET ${${PACKAGE_NAME}_LIBRARIES})
  ENDIF()

  IF (SET_INCLUDE_LIBRARY_DIRS_FROM_UPSTREAM)
     IF (FULL_INCLUDE_DIRS_SET)
       LIST(REMOVE_DUPLICATES FULL_INCLUDE_DIRS_SET)
     ENDIF()
     IF (FULL_LIBRARY_DIRS_SET)
       LIST(REMOVE_DUPLICATES FULL_LIBRARY_DIRS_SET)
     ENDIF()
  ENDIF()

  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    MESSAGE("-- " "*** Final sets of packages, libs, include dirs, and lib dirs:")
    PRINT_VAR(FULL_PACKAGE_SET)
    PRINT_VAR(FULL_LIBRARY_SET)
    PRINT_VAR(FULL_INCLUDE_DIRS_SET)
    PRINT_VAR(FULL_LIBRARY_DIRS_SET)
  ENDIF()

  #
  # C) Get the set of TPLs for this package that are enabled
  #

  # C.1) Get the set of enabled TPLs

  SET(FULL_TPL_SET "")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    LIST(APPEND FULL_TPL_SET ${${TRIBITS_PACKAGE}_LIB_REQUIRED_DEP_TPLS})
    SET(OPTIONAL_TPLS ${${TRIBITS_PACKAGE}_LIB_OPTIONAL_DEP_TPLS})
    FOREACH(TPL ${OPTIONAL_TPLS})
      # Only add if support for the optional TPL is enabled in this
      # package.  Don't just check if the TPL is enabled!
      IF(${TRIBITS_PACKAGE}_ENABLE_${TPL})
        LIST(APPEND FULL_TPL_SET ${TPL})
      ENDIF()
    ENDFOREACH()
  ENDFOREACH()
  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    PRINT_VAR(FULL_TPL_SET)
  ENDIF()

  # C.2) Sort the TPLs according to the master TPL list

  #We will use the complete list of supported tpls for the project
  #to help us create a properly ordered list of tpls.
  IF (FULL_TPL_SET)
    SET(ORDERED_FULL_TPL_SET ${FULL_TPL_SET}) 
    TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST("${${PROJECT_NAME}_REVERSE_TPLS}"
      ORDERED_FULL_TPL_SET)
  ENDIF()

  IF (TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES_DEBUG_DUMP)
    PRINT_VAR(ORDERED_FULL_TPL_SET)
  ENDIF()

  #
  # D) Get the libraries, library dirs, and the include dirs for the
  # upstream enabled TPLs
  #

  SET(${PACKAGE_NAME}_TPL_LIBRARIES "")
  SET(${PACKAGE_NAME}_TPL_INCLUDE_DIRS "")
  SET(${PACKAGE_NAME}_TPL_LIBRARY_DIRS "")
  FOREACH(TPL ${ORDERED_FULL_TPL_SET})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARIES ${TPL_${TPL}_LIBRARIES})
    LIST(APPEND ${PACKAGE_NAME}_TPL_INCLUDE_DIRS ${TPL_${TPL}_INCLUDE_DIRS})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARY_DIRS ${TPL_${TPL}_LIBRARY_DIRS})
  ENDFOREACH()

  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  #
  # E) Deal with the library rpath issues with shared libs
  #

  # Write the specification of the rpath if necessary. This is only needed if
  # we're building shared libraries.

  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${FULL_LIBRARY_DIRS_SET}")
    SET(SHARED_LIB_RPATH_COMMAND
      ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  #
  # F) Create the contents of the cmake Config.cmake file for the build tree
  #
  # Creating this file in the base dir of the package since it is possible
  # that the cmake returned path where the file was found would be useful for
  # a package, and having to dig through the hiding that is done wouldn't be
  # nice.
  #

  IF (PARSE_WRITE_CMAKE_CONFIG_FILE)

    # Custom code in configuration file.
    SET(PACKAGE_CONFIG_CODE "")
  
    # Import build tree targets into applications.
    IF(FULL_LIBRARY_SET)
      SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
"
      )
    ENDIF()

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in
      "${PARSE_WRITE_CMAKE_CONFIG_FILE}"
      )
  ENDIF()

  #
  # G) Create the export makefile for the build tree
  # 
  # This is the equivalent of the cmake version only slightly changed so that
  # it can be directly imported into a Makefile.
  #

  IF(PARSE_WRITE_EXPORT_MAKLEFILE)

    TRIBITS_LIST_TO_STRING("${FULL_LIBRARY_SET}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_FULL_LIBRARY_SET)
    TRIBITS_LIST_TO_STRING("${FULL_LIBRARY_DIRS_SET}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${FULL_INCLUDE_DIRS_SET}" "-I" MAKEFILE_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PACKAGE_NAME}_TPL_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    TRIBITS_LIBRARY_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARIES)

    TRIBITS_LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")

    TRIBITS_LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    TRIBITS_LIST_TO_STRING("${ORDERED_FULL_TPL_SET}" "" MAKEFILE_ORDERED_FULL_TPL_SET)

    # create an upper case name of the package so that we can make deprecated
    # versions of them to help people transistioning from the autotools
    # version diagnose any missed variables.
    STRING(TOUPPER ${EXPORT_FILE_VAR_PREFIX} EXPORT_FILE_VAR_PREFIX_UPPER)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in
      "${PARSE_WRITE_EXPORT_MAKLEFILE}"
      )
  ENDIF()

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
  SET(FULL_LIBRARY_DIRS_SET "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_LIB_DIR}")
  SET(FULL_INCLUDE_DIRS_SET "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")

  # Custom code in configuration file.
  SET(PACKAGE_CONFIG_CODE "")

  # Import install tree targets into applications.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"\${CMAKE_CURRENT_LIST_DIR}/../${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  IF (PARSE_WRITE_INSTALL_CMAKE_CONFIG_FILE)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
      )

  ENDIF()

  #
  # I) Write the export makefile for the install tree
  #

  IF (PARSE_WRITE_INSTALL_EXPORT_MAKLEFILE)

    # Generated Make imports must use CMAKE_INSTALL_PREFIX, rather
    # than the more platform friendly method of locating the libraries
    # and includes using the config file path above. The underlying
    # assumption here is that a generator that uses
    # CMAKE_INSTALL_PREFIX is being used.
    SET(FULL_LIBRARY_DIRS_SET ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
    SET(FULL_INCLUDE_DIRS_SET ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

    TRIBITS_LIST_TO_STRING("${FULL_LIBRARY_DIRS_SET}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${FULL_INCLUDE_DIRS_SET}" "-I" MAKEFILE_INCLUDE_DIRS)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      )
  ENDIF()

ENDFUNCTION()


#
# Set the install targets for the package config and export makefiles.
#
# The INSTALL() commands must be in a different subroutine or CMake will not
# allow you to call the rountine, even if you if() it out!
# 

FUNCTION(TRIBITS_WRITE_PROJECT_CLIENT_EXPORT_FILES_INSTALL_TARGETS PACKAGE_NAME)

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
      RENAME ${PACKAGE_NAME}Config.cmake
      )

  ENDIF()

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PACKAGE_NAME}
      )
  ENDIF()


ENDFUNCTION()


#
# Generate the ${PACKAGE_NAME}Config.cmake and/or the Makefile.export.${PACKAGE_NAME}
# for package PACKAGE_NAME.
#
# ToDo: Finish documentation!
#

FUNCTION(TRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES PACKAGE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nTRIBITS_WRITE_PACKAGE_CLIENT_EXPORT_FILES: ${PACKAGE_NAME}")
  ENDIF()

  SET(EXPORT_FILES_ARGS PACKAGE_NAME ${PACKAGE_NAME})

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    SET(WRITE_CMAKE_CONFIG_FILE
      ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake)
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("For package ${PACKAGE_NAME} creating ${WRITE_CMAKE_CONFIG_FILE}")
    ENDIF()
    APPEND_SET(EXPORT_FILES_ARGS
      WRITE_CMAKE_CONFIG_FILE "${WRITE_CMAKE_CONFIG_FILE}"
      WRITE_INSTALL_CMAKE_CONFIG_FILE)
  ENDIF()

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    SET(WRITE_EXPORT_MAKLEFILE
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PACKAGE_NAME})
    IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE("For package ${PACKAGE_NAME} creating ${WRITE_EXPORT_MAKLEFILE}")
    ENDIF()
    APPEND_SET(EXPORT_FILES_ARGS
      WRITE_EXPORT_MAKLEFILE "${WRITE_EXPORT_MAKLEFILE}"
      WRITE_INSTALL_EXPORT_MAKLEFILE)
  ENDIF()

  TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(${EXPORT_FILES_ARGS})

  TRIBITS_WRITE_PROJECT_CLIENT_EXPORT_FILES_INSTALL_TARGETS(${PACKAGE_NAME})

ENDFUNCTION()


#
# Write the outer TriBITS project configure and/or export makefiles
#
# ToDo: Finish documentation!
#

FUNCTION(TRIBITS_WRITE_PROJECT_CLIENT_EXPORT_FILES)

  TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET()

  # Reversing the package list so that libraries will be produced in order of
  # most dependent to least dependent.
  SET(PACKAGE_LIST ${${PROJECT_NAME}_SE_PACKAGES})
  LIST(REVERSE PACKAGE_LIST)


  # Loop over all packages to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_PACKAGE_SET "")
  SET(FULL_LIBRARY_SET "")
  SET(FULL_INCLUDE_DIRS_SET "")
  SET(FULL_LIBRARY_DIRS_SET "")
  FOREACH(TRIBITS_PACKAGE ${PACKAGE_LIST})
    IF(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      LIST(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
      LIST(APPEND FULL_INCLUDE_DIRS_SET ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
      LIST(APPEND FULL_LIBRARY_DIRS_SET ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()

  SET(${PROJECT_NAME}_CONFIG_LIBRARIES ${FULL_LIBRARY_SET})

  # Reversing the tpl list so that the list of tpls will be produced in
  # order of most dependent to least dependent.
  SET(TPL_LIST ${${PROJECT_NAME}_TPLS})
  LIST(REVERSE TPL_LIST)

  # Loop over all TPLs to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_TPL_SET "")
  SET(FULL_TPL_LIBRARY_SET "")
  SET(FULL_TPL_INCLUDE_DIRS_SET "")
  SET(FULL_TPL_LIBRARY_DIRS_SET "")
  FOREACH(TPL ${TPL_LIST})
    IF(TPL_ENABLE_${TPL})
      LIST(APPEND FULL_TPL_SET ${TPL})
      LIST(APPEND FULL_TPL_LIBRARY_SET ${TPL_${TPL}_LIBRARIES})
      LIST(APPEND FULL_TPL_INCLUDE_DIRS_SET ${TPL_${TPL}_INCLUDE_DIRS})
      LIST(APPEND FULL_TPL_LIBRARY_DIRS_SET ${TPL_${TPL}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()

  # it is possible that tpls are in the same directory, to keep from
  # having a very long include path or library path we will strip out
  # any duplicates. This shouldn't affect which include or library is
  # found since the first instance of any path will be the one that is
  # kept.
  LIST(REMOVE_DUPLICATES FULL_TPL_INCLUDE_DIRS_SET)
  LIST(REMOVE_DUPLICATES FULL_TPL_LIBRARY_DIRS_SET)

  SET(${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS ${FULL_TPL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS ${FULL_TPL_LIBRARY_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARIES ${FULL_TPL_LIBRARY_SET})

  #
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")

  #Config file for setting variables and finding include/library paths from the build directory
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the build tree.
  IF(FULL_LIBRARY_SET)
    LIST(SORT FULL_LIBRARY_SET)
    LIST(REMOVE_DUPLICATES FULL_LIBRARY_SET)
    EXPORT(TARGETS ${FULL_LIBRARY_SET} FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"${${TRIBITS_PACKAGE}_BINARY_DIR}/${TRIBITS_PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake )
  ENDIF()

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARIES)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    TRIBITS_LIBRARY_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARIES)

    TRIBITS_LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")

    TRIBITS_LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    TRIBITS_LIST_TO_STRING("${FULL_TPL_SET}" "" MAKEFILE_FULL_TPL_SET)

    IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
      CONFIGURE_FILE(
        ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
        ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME})
    ENDIF()
  ENDIF()

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
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_LIB_DIR}")

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries.
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the install tree.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    INSTALL(
      EXPORT ${PROJECT_NAME}
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      FILE ${PROJECT_NAME}Targets.cmake
      )
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"\${CMAKE_CURRENT_LIST_DIR}/../${TRIBITS_PACKAGE}/${TRIBITS_PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  IF (${PROJECT_NAME}_ENABLE_INSTALL_CMAKE_CONFIG_FILES)
    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      RENAME ${PROJECT_NAME}Config.cmake
      )
  ENDIF()

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
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
    SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
    SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})

    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PROJECT_NAME}
      )
  ENDIF()

  #
  # Configure the version file for ${PROJECT_NAME}
  #

  CONFIGURE_FILE(
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsConfigVersionTemplate.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake @ONLY)

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
  )

ENDFUNCTION()

