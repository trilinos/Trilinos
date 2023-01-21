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

include(TribitsGetPackageEnableStatus)


# @FUNCTION: tribits_get_sublist_enabled()
#
# Get sub-list of enabled packages
#
# Usage::
#
#   tribits_get_sublist_enabled( <enableListName>
#     <enabledSublistNameOut>  <numEnabledVarOut>)
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
#     <nondisabledListNameOut>  <numNondisabledVarOut> )
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
#     <disabledSublistNameOut>  <numDisabledVarOut>)
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
#     <nonenabledListNameOut>  <numNonenabledVarOut>)
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
