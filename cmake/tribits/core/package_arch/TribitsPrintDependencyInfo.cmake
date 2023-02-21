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


include(PrintNonemptyVarWithSpaces)


# @FUNCTION: tribits_print_initial_dependency_info()
#
# Usage::
#
#   tribits_print_initial_dependency_info()
#
# Function that prints whatever initial dependency information that is
# available that is requested by the user after the initial construction of
# the package dependency graph but **before** the call of
# `tribits_adjust_package_enables()`_.
#
# Calls:
#
# * `tribits_print_tentatively_enabled_tpls()`_
# * `tribits_dump_package_dependencies_info()`_
#
function(tribits_print_initial_dependency_info)
  tribits_print_tentatively_enabled_tpls()
  tribits_dump_package_dependencies_info()
endfunction()


# @FUNCTION: tribits_print_tentatively_enabled_tpls()
#
# Usage::
#
#   tribits_print_tentatively_enabled_tpls()
#
# Function that print the set of tentatively enabled TPLs.
#
# Does **not** modify any state!
#
function(tribits_print_tentatively_enabled_tpls)
  foreach(TPL ${${PROJECT_NAME}_DEFINED_TPLS})
    if (TPL_TENTATIVE_ENABLE_${TPL})
      message("-- Tentatively enabling TPL '${TPL}'")
      #print_var(TPL_ENABLE_${TPL})
    endif()
  endforeach()
endfunction()


# @FUNCTION: tribits_dump_package_dependencies_info()
#
# Usage:
#
#  tribits_dump_package_dependencies_info()
#
# Function that dumps (prints to STDOUT) the package dependency info if
# ``${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES==TRUE``.
#
# This function does **not** modify any state!
#
function(tribits_dump_package_dependencies_info)

  advanced_option(${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES
    "Dump the package dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  advanced_option(${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES
    "Dump the package forward dependency information."
    "${${PROJECT_NAME}_VERBOSE_CONFIGURE}" )

  message("\nPackage dependencies information:")

  tribits_print_project_list_var_and_num(DEFINED_TPLS)
  tribits_print_project_list_var_and_num(DEFINED_INTERNAL_TOPLEVEL_PACKAGES)
  tribits_print_project_list_var_and_num(DEFINED_TOPLEVEL_PACKAGES)
  tribits_print_project_list_var_and_num(DEFINED_INTERNAL_PACKAGES)
  tribits_print_project_list_var_and_num(DEFINED_PACKAGES)

  if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_PACKAGES})
      tribits_print_package_dependencies(${TRIBITS_PACKAGE})
      message("")
    endforeach()
  endif()

endfunction()


# Optionally print out a global project list var with spaces and
# unconditionally print out the number of items in that list
#
function(tribits_print_project_list_var_and_num  listVarSuffix)
  message("")
  if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    print_nonempty_var_with_spaces(${PROJECT_NAME}_${listVarSuffix}  wasPrinted)
  endif()
  print_var(${PROJECT_NAME}_NUM_${listVarSuffix})
endfunction()


# Print out dependencies for a package
#
function(tribits_print_package_dependencies  packageName)

  set(printedVar "")

  # Print deps vars

  if (printedVar)
    message("")
  endif()

  tribits_print_nonempty_package_deps_list(${packageName}  LIB  DEFINED  printedVar)
  tribits_print_nonempty_package_deps_list(${packageName}  TEST  DEFINED  printedVar)

  if (${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES)
    tribits_print_nonempty_package_forward_defined_deps_list(${packageName}  LIB
      printedVar)
    tribits_print_nonempty_package_forward_defined_deps_list(${packageName}  TEST
      printedVar)
  endif()

  if (NOT printedVar)
    message("-- ${packageName}: No defined dependencies!")
  endif()

endfunction()


# Print out a dependency list including if each item is optional or required
# ('R' or 'O')
#
# Usage::
#
#   tribits_print_nonempty_package_deps_list(<packageName> <libOrTest>
#     <printedListOut>)
#
# which prints out the list::
#
#   <packageName>_<libOrTest>_DEFINED_DEPENDENCIES
#
# if it is non-empty.
#
function(tribits_print_nonempty_package_deps_list  packageName
    libOrTest  definedOrEnabled  printedListOut
  )
  set(depsListName ${packageName}_${libOrTest}_${definedOrEnabled}_DEPENDENCIES)
  if (NOT "${${depsListName}}" STREQUAL "")
    set(lineStr "-- ${depsListName}:")
    foreach (depPkg IN LISTS ${depsListName})
      string(APPEND lineStr " ${depPkg}")
      if (${packageName}_${libOrTest}_DEP_REQUIRED_${depPkg})
        string(APPEND lineStr "[R]")
      else()
        string(APPEND lineStr "[O]")
      endif()
    endforeach()
    message("${lineStr}")
    set(${printedListOut}  TRUE  PARENT_SCOPE)
  endif()
endfunction()


# Print out a forward dependency list including if each item is optional or
# required ('R' or 'O')
#
# Usage::
#
#   tribits_print_nonempty_package_forward_defined_deps_list(<packageName> <libOrTest>
#     <printedListOut>)
#
# which prints out the list::
#
#   <packageName>_FORWARD_<libOrTest>_DEFINED_DEPENDENCIES
#
# if it is non-empty.
#
function(tribits_print_nonempty_package_forward_defined_deps_list  packageName
    libOrTest  printedListOut
  )
  set(fwdDepsListName ${packageName}_FORWARD_${libOrTest}_DEFINED_DEPENDENCIES)
  if (NOT "${${fwdDepsListName}}" STREQUAL "")
    set(lineStr "-- ${fwdDepsListName}:")
    foreach (fwdDepPkg IN LISTS ${fwdDepsListName})
      string(APPEND lineStr " ${fwdDepPkg}")
      if (${fwdDepPkg}_${libOrTest}_DEP_REQUIRED_${packageName})
        string(APPEND lineStr "[R]")
       else()
        string(APPEND lineStr "[O]")
      endif()
    endforeach()
    message("${lineStr}")
    set(${printedListOut}  TRUE  PARENT_SCOPE)
  endif()
endfunction()


#  LocalWords:  TRIBITS
