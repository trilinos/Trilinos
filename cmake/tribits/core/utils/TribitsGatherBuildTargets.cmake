# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
