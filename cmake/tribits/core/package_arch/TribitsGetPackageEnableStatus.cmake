# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# @FUNCTION: tribits_get_package_enable_status()
#
# Function that determines if a given external or internal package's enable
# status (e.g. 'ON' or 'OFF')
#
# Usage::
#
#   tribits_get_package_enable_status(<packageName>  <packageEnableOut>
#     <packageEnableVarNameOut>)
#
# On return, if non-empty, the variable ``<packageEnableOut>`` will contain
# the actual value of ``${${PROJECT_NAME}_ENABLE_<packageName>}`` or
# ``${TPL_ENABLE_<packageName>}`` or will return empty "".  If
# ``${packageName}_PACKAGE_BUILD_STATUS == "INTERNAL", then only the value of
# ``${PROJECT_NAME}_ENABLE_<packageName>`` will be considered.
#
# On return, if non-empty, the variable ``<packageEnableVarNameOut>`` will be
# either ``${${PROJECT_NAME}_ENABLE_<packageName>}`` or
# ``${TPL_ENABLE_<packageName>}``, depending on which one is used to obtain
# the value ``<packageEnableOut>``.
#
# This works for both external packages/TPLs and internal packages.
#
function(tribits_get_package_enable_status  packageName  packageEnableOut
    packageEnableVarNameOut
  )
  tribits_get_package_enable_status_assert_args("${packageName}"
    "${packageEnableOut}" "${packageEnableVarNameOut}")
  tribits_assert_package_enable_status(${packageName})
  # Determine which variable, if any, to extract enable status from
  set(packageEnableVarName "")
  if (NOT "${${PROJECT_NAME}_ENABLE_${packageName}}" STREQUAL "")
    set(packageEnableVarName  ${PROJECT_NAME}_ENABLE_${packageName})
  elseif (NOT "${TPL_ENABLE_${packageName}}" STREQUAL "")
    set(packageEnableVarName  TPL_ENABLE_${packageName})
  elseif (${packageName}_PACKAGE_BUILD_STATUS  STREQUAL  "INTERNAL")
    set(packageEnableVarName  ${PROJECT_NAME}_ENABLE_${packageName})
  elseif (${packageName}_PACKAGE_BUILD_STATUS  STREQUAL  "EXTERNAL")
    set(packageEnableVarName  TPL_ENABLE_${packageName})
  else()
    message(FATAL_ERROR "Could not determine status of package ${packageName}")
  endif()
  # Extract the enable status, set output args
  set(packageEnable ${${packageEnableVarName}})
  if (packageEnableOut)
    set(${packageEnableOut} ${packageEnable} PARENT_SCOPE)
  endif()
  if (packageEnableVarNameOut)
    set(${packageEnableVarNameOut} ${packageEnableVarName} PARENT_SCOPE)
  endif()
endfunction()


# @FUNCTION: tribits_assert_package_enable_status()
#
# Function that asserts that if both ``${PROJECT_NAME}_ENABLE_${packageName}``
# and ``TPL_ENABLE_${packageName}`` are both set to non-empty, then they must
# be the same value or this is an error.
#
# Usage::
#
#   tribits_assert_package_enable_status(<packageName>)
#
function(tribits_assert_package_enable_status  packageName)
  if ( (NOT "${${PROJECT_NAME}_ENABLE_${packageName}}" STREQUAL "")
      AND (NOT "${TPL_ENABLE_${packageName}}" STREQUAL "")
      AND (NOT "${${PROJECT_NAME}_ENABLE_${packageName}}" STREQUAL
        "${TPL_ENABLE_${packageName}}")
      )
  message(SEND_ERROR "Error, ${PROJECT_NAME}_ENABLE_${packageName}="
    "'${${PROJECT_NAME}_ENABLE_${packageName}}' !="
    " TPL_ENABLE_${packageName} = '${TPL_ENABLE_${packageName}}'")
  endif()
endfunction()
# ToDo: Create a cache var for the mode of message() above in case the user
# wants to disable the warning.


function(tribits_get_package_enable_status_assert_args  packageName  packageEnableOut
    packageEnableVarNameOut
  )
  if ("${packageName}" STREQUAL "")
    message(FATAL_ERROR "Error, packageName='' is not allowed!")
  endif()
  if ("${packageEnableOut}" STREQUAL "" AND "${packageEnableVarNameOut}" STREQUAL "")
    message(FATAL_ERROR "Error, both packageEnableOut='' and"
      " packageEnableVarNameOut='' is not allowed!")
  endif()
endfunction()
