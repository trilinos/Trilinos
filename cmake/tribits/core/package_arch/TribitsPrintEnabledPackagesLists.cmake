# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsGetPackageSublists)


# @FUNCTION: tribits_print_enables_before_adjust_package_enables()
#
# Call this to print the package enables before calling
# tribits_adjust_package_enables().
#
# Usage::
#
#   tribits_print_enables_before_adjust_package_enables()
#
function(tribits_print_enables_before_adjust_package_enables)
  tribits_print_toplevel_package_list_enable_status(
    "\nExplicitly enabled top-level packages on input (by user)"
    INTERNAL ON NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nExplicitly enabled packages on input (by user)" INTERNAL ON NONEMPTY)
  tribits_print_toplevel_package_list_enable_status(
    "\nExplicitly disabled top-level packages on input (by user or by default)"
    INTERNAL OFF NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nExplicitly disabled packages on input (by user or by default)"
    INTERNAL OFF NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nExplicitly enabled external packages/TPLs on input (by user)" EXTERNAL ON NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nExplicitly disabled external packages/TPLs on input (by user or by default)"
    EXTERNAL OFF NONEMPTY)
  tribits_print_package_build_status("\nInitial package build status:\n"
    "-- Initial: " PRINT_ALL)
endfunction()


# @FUNCTION: tribits_print_enables_after_adjust_package_enables()
#
# Call this to print the package enables before calling
# tribits_adjust_package_enables().
#
# Usage::
#
#   tribits_print_enables_after_adjust_package_enables()
#
function(tribits_print_enables_after_adjust_package_enables)
  tribits_print_toplevel_package_list_enable_status(
   "\nFinal set of enabled top-level packages" INTERNAL ON NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nFinal set of enabled packages" INTERNAL ON NONEMPTY)
  tribits_print_toplevel_package_list_enable_status(
    "\nFinal set of non-enabled top-level packages" INTERNAL OFF INCLUDE_EMPTY)
  tribits_print_package_list_enable_status(
    "\nFinal set of non-enabled packages" INTERNAL OFF INCLUDE_EMPTY)
  tribits_print_toplevel_package_list_enable_status(
    "\nFinal set of enabled top-level external packages/TPLs" EXTERNAL ON NONEMPTY)
  tribits_print_package_list_enable_status(
    "\nFinal set of enabled external packages/TPLs" EXTERNAL ON NONEMPTY)
  tribits_print_toplevel_package_list_enable_status(
    "\nFinal set of non-enabled top-level external packages/TPLs" EXTERNAL OFF INCLUDE_EMPTY)
  tribits_print_package_list_enable_status(
    "\nFinal set of non-enabled external packages/TPLs" EXTERNAL OFF INCLUDE_EMPTY)
  tribits_print_package_build_status("\nFinal package build status (enabled only):\n"
    "-- Final: " PRINT_ONLY_ENABLED)
endfunction()


# Function that prints the current set of enabled internal or external
# top-level packages
#
function(tribits_print_toplevel_package_list_enable_status
    docString  internalOrExternal  enabledFlag  enableEmptyStatus
  )
  tribits_print_packages_list_enable_status_from_var(
    ${PROJECT_NAME}_DEFINED_TOPLEVEL_PACKAGES
    "${docString}" "${internalOrExternal}" ${enabledFlag} ${enableEmptyStatus} )
endfunction()


# Prints the current set of enabled/disabled internal packages
#
function(tribits_print_package_list_enable_status
    docString  internalOrExternal  enabledFlag  enableEmptyStatus
  )
  tribits_print_packages_list_enable_status_from_var(
    ${PROJECT_NAME}_DEFINED_PACKAGES
    "${docString}" "${internalOrExternal}" ${enabledFlag} ${enableEmptyStatus} )
endfunction()


# Print the current set of enabled/disabled packages given input list of
# packages
#
function(tribits_print_packages_list_enable_status_from_var  packageListVarName
    docString  internalOrExternal  enabledFlag  enableEmptyStatus
  )
  tribits_filter_package_list_from_var(${packageListVarName}
    "${internalOrExternal}"  "${enabledFlag}"  ${enableEmptyStatus} packageSublist)
  tribits_print_prefix_string_and_list("${docString}"  packageSublist)
endfunction()


# Print out a list with white-space separators with an initial doc string
#
function(tribits_print_prefix_string_and_list  docString   listNameToPrint)
  string(REPLACE ";" " " listToPrintStr "${${listNameToPrint}}")
  list(LENGTH  ${listNameToPrint}  numElements)
  if (numElements GREATER "0")
    message("${docString}:  ${listToPrintStr} ${numElements}")
  else()
    message("${docString}:  ${numElements}")
  endif()
endfunction()


# Print out the `<Package>_PACKAGE_BUILD_STATUS` vars
#
function(tribits_print_package_build_status  summaryHeaderStr  prefixStr
    printEnabledOpt
  )
  if (printEnabledOpt STREQUAL "PRINT_ALL")
    set(printAll ON)
  elseif (printEnabledOpt STREQUAL "PRINT_ONLY_ENABLED")
    set(printAll OFF)
  else()
    message(FATAL_ERROR
      "Error, invalid value for printEnabledOpt='${printEnabledOpt}'")
  endif()

  if (${PROJECT_NAME}_DUMP_PACKAGE_BUILD_STATUS)
    message("${summaryHeaderStr}")
    foreach(packageName  IN LISTS  ${PROJECT_NAME}_DEFINED_PACKAGES)
    tribits_get_package_enable_status(${packageName}  packageEnable  "")
      if (printAll OR packageEnable)
        message("${prefixStr}${packageName}_PACKAGE_BUILD_STATUS="
	  "${${packageName}_PACKAGE_BUILD_STATUS}")
      endif()
    endforeach()
  endif()
endfunction()
