# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @MACRO: tribits_pkg_export_cache_var()
#
# Macro that registers a package-level cache var to be exported in the
# ``<Package>Config.cmake`` file
#
# Usage::
#
#   tribits_pkg_export_cache_var(<cacheVarName>)
#
# where ``<cacheVarName>`` must be the name of a cache variable (or an error
# will occur).
#
# NOTE: This will also export this variable to the
# ``<Package><Spkg>Config.cmake`` file for every enabled subpackage (if this
# is called from a ``CMakeLists.txt`` file of a top-level package that has
# subpackages).  That way, any top-level package cache vars are provided by
# any of the subpackages' ``<Package><Spkg>Config.cmake`` files.
#
macro(tribits_pkg_export_cache_var  cacheVarName)
  if (DEFINED ${PACKAGE_NAME}_PKG_VARS_TO_EXPORT)
    # Assert this is a cache var
    get_property(cacheVarIsCacheVar  CACHE ${cacheVarName} PROPERTY  VALUE  SET)
    if (NOT cacheVarIsCacheVar)
      message(SEND_ERROR
        "ERROR: The variable ${cacheVarName} is NOT a cache var and cannot"
        " be exported!")
    endif()
    # Add to the list of package cache vars to export
    append_global_set(${PACKAGE_NAME}_PKG_VARS_TO_EXPORT
      ${cacheVarName})
  endif()
endmacro()


# @MACRO: tribits_assert_cache_and_local_vars_same_value()
#
# Asset that a cache variable and a possible local variable (if it exists)
# have the same value.
#
# Usage::
#
#   tribits_assert_cache_and_local_vars_same_value(<cacheVarName>)
#
# If the local var ``<cacheVarName>`` and the cache var ``<cacheVarName>``
# both exist but have different values, then ``message(SEND_ERROR ...)`` is
# called with an informative error message.
#
macro(tribits_assert_cache_and_local_vars_same_value  cacheVarName)
  set(cacheVarValue "$CACHE{${cacheVarName}}")
  set(localValue "${${cacheVarName}}")
  if (NOT localValue STREQUAL cacheVarValue)
    message_wrapper(SEND_ERROR "ERROR: The cache variable ${cacheVarName} with the"
      " cache var value '${cacheVarValue}' is not the same value as the local"
      " variable ${cacheVarName} with value '${localValue}'!")
  endif()
endmacro()


# Function that sets up data-structures for package-level cache var to be
# exported
#
function(tribits_pkg_init_exported_vars  PACKAGE_NAME_IN)
  global_set(${PACKAGE_NAME_IN}_PKG_VARS_TO_EXPORT "")
endfunction()


# Function that injects set() statements for a package's exported cache vars into
# a string.
#
# This is used to create set() statements to be injected into a package's
# ``<Package>Config.cmake`` file.
#
function(tribits_pkg_append_set_commands_for_exported_vars  packageName
    configFileStrInOut
  )
  set(configFileStr "${${configFileStrInOut}}")
  if (NOT "${${packageName}_PARENT_PACKAGE}" STREQUAL "")
    foreach(exportedCacheVar IN LISTS ${${packageName}_PARENT_PACKAGE}_PKG_VARS_TO_EXPORT)
      tribits_assert_cache_and_local_vars_same_value(${exportedCacheVar})
      string(APPEND configFileStr
        "set(${exportedCacheVar} \"${${exportedCacheVar}}\")\n")
    endforeach()
  endif()
  foreach(exportedCacheVar IN LISTS ${packageName}_PKG_VARS_TO_EXPORT)
    #tribits_assert_cache_and_local_vars_same_value(${exportedCacheVar})
    string(APPEND configFileStr
      "set(${exportedCacheVar} \"${${exportedCacheVar}}\")\n")
  endforeach()
  set(${configFileStrInOut} "${configFileStr}" PARENT_SCOPE)
endfunction()
