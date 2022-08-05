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


# Function to get a list all of buildsystem targets of the type matching the
# given list starting in a subdir and in below subdirs.
#
# Usage::
#
#   tribits_get_all_build_targets_including_in_subdirs(
#     <subdir> <targetTypesList> <targetsListVarOut> )
#
# On output the variable ``<targetsListVarOut>`` will contain a list of all of
# the build targets match match the targets types listed in
# ``<targetTypesList>`` which is a quoted list of the form
# ``"<targettype0>;<targettype2>;..."``.  For example, to get all of the
# real library and executable targets, run::
#
#   tribits_get_all_build_targets_including_in_subdirs(
#     <subdir>
#     "STATIC_LIBRARY;SHARED_LIBRARY;EXECUTABLE"
#     <targetsListVarOut> )
#
function(tribits_get_all_build_targets_including_in_subdirs  srcdir  targetTypesList
    targetsListVarOut
  )

  set(targetsList "")

  # Recurse into subdirectories.
  get_property(dirs DIRECTORY "${srcdir}" PROPERTY SUBDIRECTORIES)
  foreach(d IN LISTS dirs)
    tribits_get_all_build_targets_including_in_subdirs(${d} "${targetTypesList}"
      targetsSubdirList)
    list(APPEND targetsList ${targetsSubdirList})
  endforeach()

  # Get the targets from this directory.
  get_property(allTargetsThisDir DIRECTORY "${srcdir}" PROPERTY BUILDSYSTEM_TARGETS)
  tribits_filter_only_build_targets(allTargetsThisDir "${targetTypesList}"
    buildTargetsThisDir)
  list(APPEND targetsList ${buildTargetsThisDir})

  # Return
  set(${targetsListVarOut} ${targetsList} PARENT_SCOPE)

endfunction()


function(tribits_filter_only_build_targets targetListInVar targetTypesList
    targetListOutVar
  )

  set(targetListOut "")

  foreach (target IN LISTS ${targetListInVar})
    get_property(targetType TARGET ${target} PROPERTY TYPE)
    foreach (desiredTargetType IN ITEMS ${targetTypesList})
      if (desiredTargetType STREQUAL targetType)
        list(APPEND targetListOut ${target})
        break()
      endif()
    endforeach()
  endforeach()

  set(${targetListOutVar} ${targetListOut} PARENT_SCOPE)

endfunction()
