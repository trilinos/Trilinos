# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsGetPackageEnableStatus)


# @FUNCTION: tribits_filter_package_list_from_var()
#
# Filter a list of packages based on several criteria including
# internal/external, enable status (with empty or non-empty)
#
# Usage::
#
#  tribits_filter_package_list_from_var( <packageListVarName>
#    <internalOrExternal>  <enabledFlag>  <enableEmptyStatus>
#    <packageSublistOut> )
#
# Where:
#
# * ``<packageListVarName>``: Name of input list var of packages
# * ``<internalOrExternal>``: ``INTERNAL``, ``EXTERNAL`` or "" (i.e. both
#   INTERNAL and EXTERNAL packages) (matches
#   ``<Package>_PACKAGE_BUILD_STATUS``)
# * ``<enabledFlag>``: ``ON`` for elements that match ``TRUE`` , ``OFF``
#   for elements that match ``FALSE``
# * ``<enableEmptyStatus>``: Enable status is ``NONEMPTY`` (i.e. must have a
#   value) or ``INCLUDE_EMPTY`` (matches those with a value and with empty
#   "").
# * ``<packageSublistOut>``: The name of the var that will store the filtered
#   sublist of packages.
#
function(tribits_filter_package_list_from_var  packageListVarName
    internalOrExternal  enabledFlag  enableEmptyStatus  packageSublistOut
  )
  tribits_assert_include_empty_param(${enableEmptyStatus})
  tribits_get_sublist_internal_external(${packageListVarName} "${internalOrExternal}"
    packageListTmp "")
  if (enabledFlag  AND  (enableEmptyStatus STREQUAL "NONEMPTY"))
    tribits_get_sublist_enabled(packageListTmp
      packageSublist  "")
  elseif (enabledFlag  AND  (enableEmptyStatus STREQUAL "INCLUDE_EMPTY"))
    tribits_get_sublist_nondisabled(${packageListTmp}  packageSublist  "")
  elseif (NOT enabledFlag  AND  (enableEmptyStatus STREQUAL "NONEMPTY"))
    tribits_get_sublist_disabled(packageListTmp
      packageSublist  "")
  elseif (NOT enabledFlag  AND  (enableEmptyStatus STREQUAL "INCLUDE_EMPTY"))
    tribits_get_sublist_nonenabled(packageListTmp
      packageSublist  "")
  else()
    message(FATAL_ERROR "Should never get here!")
  endif()
  set(${packageSublistOut} ${packageSublist} PARENT_SCOPE)
endfunction()


# @FUNCTION: tribits_get_sublist_enabled()
#
# Get sub-list of enabled packages
#
# Usage::
#
#   tribits_get_sublist_enabled( <enableListName>
#     <enabledSublistNameOut>  [<numEnabledVarOut>] )
#
# On output, ``<enabledSublistNameOut>`` contains the sublist of entries in
# ``<enableListName>`` which evaluate to ``TRUE`` in an ``if ()`` statement.
#
function(tribits_get_sublist_enabled  enableListName
    enabledSublistNameOut  numEnabledVarOut
  )
  set(enabledSublist)
  foreach(pkgName  IN LISTS  ${enableListName})
    tribits_get_package_enable_status(${pkgName}  enableVal  "")
    if (enableVal)
      list(APPEND  enabledSublist  ${pkgName})
    endif()
  endforeach()
  list(LENGTH  enabledSublist  numEnabled)
  set(${enabledSublistNameOut} ${enabledSublist} PARENT_SCOPE)
  if (numEnabledVarOut)
    set(${numEnabledVarOut} ${numEnabled} PARENT_SCOPE)
  endif()
endfunction()


# @FUNCTION: tribits_get_sublist_nondisabled()
#
# Get sub-list of non-disabled packages
#
# Usage::
#
#   tribits_get_sublist_nondisabled( <enableListName>
#     <nondisabledListNameOut>  [<numNondisabledVarOut>] )
#
# On output, ``<nondisabledListNameOut>`` contains the sublist of entries from
# ``<enableListName>`` for which evaluate to ``TRUE`` or empty ``""`` in an
# ``if ()`` statement.
#
function(tribits_get_sublist_nondisabled  enableListName
    nondisabledListNameOut  numNondisabledVarOut
  )
  set(nondisabledList "")
  foreach(pkgName  IN LISTS  ${enableListName})
    tribits_get_package_enable_status(${pkgName}  enableVal  "")
    if (enableVal OR "${enableVal}" STREQUAL "")
      list(APPEND  nondisabledList  ${pkgName})
    endif()
  endforeach()
  list(LENGTH  nondisabledList  numNondisabled)
  set(${nondisabledListNameOut} ${nondisabledList} PARENT_SCOPE)
  if (numNondisabledVarOut)
    set(${numNondisabledVarOut} ${numNondisabled} PARENT_SCOPE)
  endif()
endfunction()


# @FUNCTION: tribits_get_sublist_disabled()
#
# Get sub-list of disabled packages
#
# Usage::
#
#   tribits_get_sublist_disabled( <enableListName>
#     <disabledSublistNameOut>  [<numDisabledVarOut>] )
#
# On output, ``<disabledSublistNameOut>`` contains the sublist of entries
# ``<enableListName>`` which evaluate to ``FALSE`` and is not empty ``""`` in
# an ``if ()`` statement.
#
function(tribits_get_sublist_disabled  enableListName
    disabledSublistNameOut  numDisabledVarOut
  )
  set(disabledSublist "")
  foreach(pkgName  IN LISTS  ${enableListName})
    tribits_get_package_enable_status(${pkgName}  enableVal  "")
    if ((NOT enableVal) AND (NOT "${enableVal}" STREQUAL ""))
      list(APPEND  disabledSublist  ${pkgName})
    endif()
  endforeach()
  list(LENGTH  disabledSublist  numDisabled)
  set(${disabledSublistNameOut} ${disabledSublist} PARENT_SCOPE)
  if (numDisabledVarOut)
    set(${numDisabledVarOut} ${numDisabled} PARENT_SCOPE)
  endif()
endfunction()


# @FUNCTION: tribits_get_sublist_nonenabled()
#
# Get sub-list of non-enabled entries
#
# Usage::
#
#   tribits_get_sublist_nonenabled( <enableListName>
#     <nonenabledListNameOut>  [<numNonenabledVarOut>] )
#
# On output, ``<nonenabledListNameOut>`` contains the subset of entries in
# ``<enableListName>`` that evaluate to ``FALSE`` (which can also be empty
# ``""``) in an ``if ()`` statement.
#
function(tribits_get_sublist_nonenabled  enableListName
    nonenabledListNameOut  numNonenabledVarOut
  )
  set(nonenabledList "")
  foreach(pkgName  IN LISTS  ${enableListName})
    tribits_get_package_enable_status(${pkgName}  enableVal  "")
    if (NOT enableVal)
      list(APPEND  nonenabledList  ${pkgName})
    endif()
  endforeach()
  list(LENGTH  nonenabledList  numNonenabled)
  set(${nonenabledListNameOut} ${nonenabledList} PARENT_SCOPE)
  if (numNonenabledVarOut)
    set(${numNonenabledVarOut} ${numNonenabled} PARENT_SCOPE)
  endif()
endfunction()


# @FUNCTION: tribits_get_sublist_internal_external()
#
# Get sub-list of packages that are INTERNAL, EXTERNAL, or either.
#
# Usage::
#
#   tribits_get_sublist_internal_external( <inputPackageListName>  <internalOrExternal>
#     <sublistNameOut>  [<sizeSublistOut>] )
#
# where:
#
#   * `<internalOrExternal>` is either `INTERNAL`, `EXTERNAL` or empty "".
#
# On output, ``<sublistNameOut>`` contains the sublist of entries in
# ``<inputPackageListName>`` which are either `INTERNAL` or `EXTERNAL` or both
# (if `<internalOrExternal>` is "") based on `<Package>_PACKAGE_BUILD_STATUS`
# for each element package name.
#
function(tribits_get_sublist_internal_external  inputPackageListName  internalOrExternal
    sublistNameOut  sizeSublistOut
  )
  tribits_assert_internal_or_external_arg("${internalOrExternal}")
  if (NOT internalOrExternal STREQUAL "")
    set(sublist "")
    foreach(pkgName  IN LISTS  ${inputPackageListName})
      if (${pkgName}_PACKAGE_BUILD_STATUS STREQUAL internalOrExternal)
        list(APPEND sublist ${pkgName})
      endif()
    endforeach()
  else()
    set(sublist ${${inputPackageListName}})
  endif()
  set(${sublistNameOut} ${sublist} PARENT_SCOPE)
  if (numEnabledVarOut)
    list(LENGTH  sublist  sizeSublist)
    set(${sizeSublistOut} ${sizeSublist} PARENT_SCOPE)
  endif()
endfunction()


function(tribits_assert_internal_or_external_arg  internalOrExternal)
  set(allowedValuesList "INTERNAL" "EXTERNAL" "")
  if (NOT internalOrExternal IN_LIST allowedValuesList)
    message(FATAL_ERROR "ERROR, the argument internalOrExternal="
      "'${internalOrExternal}' is not in the list of allowed values"
      ${allowedValuesList} )
  endif()
endfunction()


function(tribits_assert_include_empty_param  enableEmptyStatus)
  if ((NOT enableEmptyStatus STREQUAL "NONEMPTY") AND
      (NOT enableEmptyStatus STREQUAL "INCLUDE_EMPTY")
    )
    message(FATAL_ERROR "Error, argument enableEmptyStatus='${enableEmptyStatus}'"
      " is invalid!  Use 'NONEMPTY' or 'INCLUDE_EMPTY'.")
  endif()
endfunction()
