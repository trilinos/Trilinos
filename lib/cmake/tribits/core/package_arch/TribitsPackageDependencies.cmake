# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


################################################################################
#
# TribitPackageDependencies.cmake
#
# Contains macros and functions related to defining internal package and
# external package/TPL dependencies with a little bit of logic.
#
################################################################################


include_guard()

include("${CMAKE_CURRENT_LIST_DIR}/../common/TribitsCMakePolicies.cmake"
  NO_POLICY_SCOPE)

include(TribitsParseArgumentsHelpers)
include(MessageWrapper)


# @MACRO: tribits_extpkg_define_dependencies()
#
# Macro called from inside of a `FindTPL<tplName>Dependencies.cmake`_ file to
# define the direct upstream dependencies an external package/TPL.
#
# Usage::
#
#   tribits_extpkg_define_dependencies( <tplName>
#     DEPENDENCIES <upstreamTpl_0> <upstreamTpl_1>:<vis_1> ... )
#
# The listed upstream dependencies ``<upstreamTpl_i>`` are other external
# packages/TPLs listed before this external packages/TPL ``<tplName>`` in a
# `<repoDir>/TPLsList.cmake`_ file.  Each upstream dependency can include a
# visibility specification ``<vis_i>`` that can be added to the dependency
# using a colon ``:`` with ``<upstreamTpl_1>:<vis_1>`` where ``<vis_i>`` can
# take the allowed values:
#
# * ``PUBLIC``: The include directories along with the libraries will be
#   propagated downstream to clients of ``<tplName>``.
#
# * ``PRIVATE``: Only the libraries will be propagated downstream to clients
#   of ``<tplName>``.
#
# If ``<vis_i>`` is not specified, then the default is ``PRIVATE``.  (If a
# package needs the include directories from some external package/TPL, then
# it should list that external package/TPL as a direct dependency and not
# expect to get include directories from indirect dependencies.)
#
macro(tribits_extpkg_define_dependencies
    tplName
  )

  # Parse arguments
  cmake_parse_arguments(
     PARSE "" "" # prefix, options, one_value_keywords
     "DEPENDENCIES"  # multi_value_keywords
     ${ARGN}
     )
  tribits_check_for_unparsed_arguments(PARSE)
  tribits_assert_parse_arg_one_or_more_values(PARSE  DEPENDENCIES)

  set(${tplName}_LIB_DEFINED_DEPENDENCIES  ${PARSE_DEPENDENCIES}  CACHE STRING
    "List of all potential dependencies for external package/TPL '${tplName}'")
  mark_as_advanced(${tplName}_LIB_DEFINED_DEPENDENCIES)

endmacro()
#
# NOTE: Above, we use a cache var for ${tplName}_LIB_DEFINED_DEPENDENCIES to allow
# the user to override what dependencies a TPL can depend on.  Since this does
# not depend on what actual TPLs are enabled, it is okay to set this as a
# cache var.  As with any genetic change to a CMakeLists.txt file, you always
# need to blow a way the cache and configure from scratch to get a correct
# configuration.  (Only some types of changes to CMakeLists.txt files allow
# for correct reconfiguration.)


# @MACRO: tribits_extpkg_setup_enabled_dependencies()
#
# Usage::
#
#   tribits_extpkg_setup_enabled_dependencies(<externalPkgName>)
#
# Macro that sets up the list of enabled external package/TPL dependencies
#
# Takes the list ``<externalPkgName>_LIB_DEFINED_DEPENDENCIES`` and sets the
# default entries of the non-cache var
# ``<externalPkgName>_LIB_ENABLED_DEPENDENCIES``.  However, if
# ``${<externalPkgName>_LIB_ENABLED_DEPENDENCIES`` is non-empty when this
# macro is called, then it will not be changed.  That allows the user to
# override the list of enabled TPL dependencies in the cache.  This also sets
# the non-cache vars ``<externalPkgName>_ENABLE_<upstsreamPkgName>=ON`` for
# each enabled package listed in
# ``<externalPkgName>_LIB_ENABLED_DEPENDENCIES`` and to ``OFF`` for each
# ``<upstsreamPkgName>`` listed in
# ``<externalPkgName>_LIB_DEFINED_DEPENDENCIES`` but not in
# ``<externalPkgName>_LIB_ENABLED_DEPENDENCIES``.
#
macro(tribits_extpkg_setup_enabled_dependencies  externalPkgName)

  set(libEnabledDependencies "")

  if (TPL_ENABLE_${externalPkgName})
    foreach(upstreamPkgEntry  IN  LISTS ${externalPkgName}_LIB_DEFINED_DEPENDENCIES)
      tribits_extpkg_get_dep_name_and_vis(${upstreamPkgEntry}
        upstreamPkgName  upstreamPkgVis)
      if (TPL_ENABLE_${upstreamPkgName})
        list(APPEND libEnabledDependencies ${upstreamPkgEntry})
      else()
        set(${externalPkgName}_ENABLE_${upstreamPkgName} OFF)
      endif()
    endforeach()
  endif()

  if ("${${externalPkgName}_LIB_ENABLED_DEPENDENCIES}" STREQUAL "")
    # Only set of not already set as a cache var, for example, by the user
    set(${externalPkgName}_LIB_ENABLED_DEPENDENCIES  ${libEnabledDependencies})
  endif()

  foreach(upstreamPkgEntry  IN  LISTS ${externalPkgName}_LIB_ENABLED_DEPENDENCIES)
    tribits_extpkg_get_dep_name_and_vis( ${upstreamPkgEntry}
      upstreamPkgName  upstreamPkgVis)
    set(${externalPkgName}_ENABLE_${upstreamPkgName} ON)
  endforeach()

endmacro()
#
# NOTE: Above we set ${externalPkgName}_LIB_ENABLED_DEPENDENCIES as a regular
# project-level variable, not as a cache var.  That is because we want it to
# update when the set of enabled TPLs changes without having to reconfigure
# from scratch.  However, we use the if statement to allow the user to
# override the default logic on the dependencies for an enabled TPL.


# @FUNCTION: tribits_extpkg_get_dep_name_and_vis()
#
# Extract ``<PkgName>`` and ``<Vis>`` from ``<PkgName>[:<Vis>]`` input with
# default ``<Vis>`` of ``PRIVATE``.
#
# Usage::
#
#   tribits_extpkg_get_dep_name_and_vis(
#     <upstreamTplDepEntry> <upstreamTplDepNameOut> <upstreamTplDepVisOut>)
#
function(tribits_extpkg_get_dep_name_and_vis
    upstreamTplDepEntry  upstreamTplDepNameOut  upstreamTplDepVisOut
  )
  # Split on ':' to get <PkgName>[:<Vis>]
  string(REPLACE ":" ";" upstreamTplAndVisList  "${upstreamTplDepEntry}")
  list(LENGTH upstreamTplAndVisList upstreamTplAndVisListLen)
  # Validate one or two entries only
  if (upstreamTplAndVisListLen GREATER 2)
    math(EXPR numColons "${upstreamTplAndVisListLen}-1")
    message_wrapper(FATAL_ERROR
      "ERROR: '${upstreamTplDepEntry}' has ${numColons} ':' but only 1 is allowed!")
  endif()
  # Set <PkgName>
  list(GET upstreamTplAndVisList 0 upstreamTplDepName)
  # Set <Vis>
  if (upstreamTplAndVisListLen EQUAL 2)
    list(GET upstreamTplAndVisList 1 upstreamTplDepVis)
  else()
    set(upstreamTplDepVis PRIVATE)
  endif()
  # Assert valid <Vis>
  set(validVisValues "PUBLIC" "PRIVATE")
  if (NOT  upstreamTplDepVis  IN_LIST  validVisValues)
    message_wrapper(FATAL_ERROR
      "ERROR: '${upstreamTplDepEntry}' has invalid visibility '${upstreamTplDepVis}'."
      "  Only 'PUBLIC' or 'PRIVATE' allowed!")
    return()  # Only executed in unit-test mode!
  endif()
  # Set outputs
  set(${upstreamTplDepNameOut} ${upstreamTplDepName} PARENT_SCOPE)
  set(${upstreamTplDepVisOut} ${upstreamTplDepVis} PARENT_SCOPE)
endfunction()
