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


include(TribitsProcessPackagesAndDirsLists)
include(TribitsAddOptionAndDefine)
include(TribitsGeneralMacros)
include(TribitsPrintEnabledPackagesLists)
include(TribitsPackageDependencies)

include(AdvancedOption)
include(AdvancedSet)
include(AppendStringVar)
include(CMakeBuildTypesList)
include(FindListElement)
include(GlobalNullSet)
include(PrintNonemptyVar)
include(PrintVarWithSpaces)
include(PrintNonemptyVarWithSpaces)
include(PrintVar)
include(RemoveGlobalDuplicates)
include(SetDefault)
include(MessageWrapper)
include(DualScopeSet)
include(CMakeParseArguments)


#
# Private helper macros
#


function(tribits_private_print_disable
  ENABLE_BEING_DISABLED_VAR_NAME  PACKAGE_WITH_SOMETHING_BEING_DISABLED
  DEP_TYPE_STR  THING_DISALBED_TYPE  THING_DISABLED_NAME
  )
  #print_var(${ENABLE_BEING_DISABLED_VAR_NAME})
  if (${ENABLE_BEING_DISABLED_VAR_NAME})
    if (${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES)
      message(
        " ***\n"
        " *** NOTE: Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
        " which was '${${ENABLE_BEING_DISABLED_VAR_NAME}}' because"
        " ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has"
        " a required ${DEP_TYPE_STR} dependence on disabled"
        " ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}"
        " but ${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON!\n"
        " ***\n"
        )
    else()
      message(FATAL_ERROR
        " ***\n"
        " *** ERROR: Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
        " which was '${${ENABLE_BEING_DISABLED_VAR_NAME}}' because"
        " ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has"
        " a required ${DEP_TYPE_STR} dependence on disabled"
        " ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}!\n"
        " ***\n"
        )
    endif()
  else()
    message("-- "
      "Setting ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
      " because ${PACKAGE_WITH_SOMETHING_BEING_DISABLED} has a required ${DEP_TYPE_STR}"
      " dependence on disabled ${THING_DISALBED_TYPE} ${THING_DISABLED_NAME}")
  endif()
endfunction()


macro(tribits_private_disable_tpl_required_package_enable
  TPL_NAME  PACKAGE_NAME  LIBRARY_DEP
  )

  #message("TRIBITS_PRIVATE_DISABLE_TPL_REQUIRED_PACKAGE_ENABLE"
  #  " ${TPL_NAME} ${PACKAGE_NAME} ${LIBRARY_DEP}")

  # Only turn off PACKAGE_NAME libraries and test/examples if it
  # is currently enabled or could be enabled.

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
     OR ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} STREQUAL ""
     )

    if ("${LIBRARY_DEP}" STREQUAL "TRUE")

      tribits_private_print_disable(
        ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} ${PACKAGE_NAME} "library"
        "TPL" ${TPL_NAME}
        )

      set(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} OFF)

    else()

      set(DEP_TYPE_STR "test/example")

      if (${PACKAGE_NAME}_ENABLE_TESTS
        OR "${${PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
        )
        tribits_private_print_disable(
          ${PACKAGE_NAME}_ENABLE_TESTS ${PACKAGE_NAME} "${DEP_TYPE_STR}"
          "TPL" ${TPL_NAME}
          )
        set(${PACKAGE_NAME}_ENABLE_TESTS OFF)
      endif()

      if (${PACKAGE_NAME}_ENABLE_EXAMPLES
        OR "${${PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        tribits_private_print_disable(
          ${PACKAGE_NAME}_ENABLE_EXAMPLES ${PACKAGE_NAME} "${DEP_TYPE_STR}"
          "TPL" ${TPL_NAME}
          )
        set(${PACKAGE_NAME}_ENABLE_EXAMPLES OFF)
      endif()

      # NOTE: We can't assert that ${PACKAGE_NAME}_ENABLE_TESTS or
      # ${PACKAGE_NAME}_ENABLE_EXAMPLES exists yet because
      # tribits_set_up_optional_package_enables_and_cache_vars() which defines
      # them is not called until after the final package enables are set.

    endif()

  endif()

endmacro()


function(tribits_private_print_disable_required_package_enable
  PACKAGE_NAME  PACKAGE_ENABLE_SOMETHING_VAR_NAME  FORWARD_DEP_PACKAGE_NAME
  DEP_TYPE_STR
  )
  tribits_private_print_disable(
    ${PACKAGE_ENABLE_SOMETHING_VAR_NAME} ${FORWARD_DEP_PACKAGE_NAME}
    "${DEP_TYPE_STR}" "package" ${PACKAGE_NAME} )
endfunction()


macro(tribits_private_disable_required_package_enables
  FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME LIBRARY_DEP
  )

  #message("TRIBITS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES"
  #  " ${FORWARD_DEP_PACKAGE_NAME} ${LIBRARY_DEP}")

  # Only turn off FORWARD_DEP_PACKAGE libraries and test/examples if it
  # is currently enabled or could be enabled

  assert_defined(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
  if (${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}
     OR ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} STREQUAL ""
     )

    if ("${LIBRARY_DEP}" STREQUAL "TRUE")

      tribits_private_print_disable_required_package_enable(
        ${PACKAGE_NAME} ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}
        ${FORWARD_DEP_PACKAGE_NAME} "library" )

      set(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} OFF)

    else()

      set(DEP_TYPE_STR "test/example")

      if (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS
        OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
        )
        tribits_private_print_disable_required_package_enable(
          ${PACKAGE_NAME} ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS
          ${FORWARD_DEP_PACKAGE_NAME} "${DEP_TYPE_STR}" )
        set(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS OFF)
      endif()

      if (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES
        OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        tribits_private_print_disable_required_package_enable(
          ${PACKAGE_NAME} ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES
          ${FORWARD_DEP_PACKAGE_NAME} "${DEP_TYPE_STR}" )
        set(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES OFF)
      endif()

    endif()

  endif()

endmacro()


macro(tribits_private_disable_optional_package_enables
  FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME
  )

  #message("TRIBITS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES"
  #  " ${FORWARD_DEP_PACKAGE_NAME} ${PACKAGE_NAME}")
  #message("-- " "${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} = ${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}")

  #assert_defined(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME})
  if (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OR "${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}" STREQUAL "")
    # Always disable the conditional enable but only print the message if the package is enabled.
    #message("--  Disable ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} ...")
    if (${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
      if (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME})  # is explicitly try already!
        message("-- "
          "NOTE: Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}=OFF"
          " which was ${${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}}"
          " because ${FORWARD_DEP_PACKAGE_NAME} has an optional library dependence"
          " on disabled package ${PACKAGE_NAME}")
      else()  # Not explicitly set
        message("-- "
          "Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}=OFF"
          " because ${FORWARD_DEP_PACKAGE_NAME} has an optional library dependence"
          " on disabled package ${PACKAGE_NAME}")
      endif()
    endif()
    set(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OFF)
  endif()

endmacro()


# Macro that disabled a packages if its required upstream TPL is disabled..
#
macro(tribits_disable_package_if_tpl_disabled  TRIBITS_PACKAGE)

  foreach(TPL_NAME ${${TRIBITS_PACKAGE}_LIB_REQUIRED_DEP_TPLS})
    if ( (NOT TPL_ENABLE_${TPL_NAME}) AND
      (NOT "${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
      )
      tribits_private_disable_tpl_required_package_enable(
        ${TPL_NAME}  ${TRIBITS_PACKAGE}  TRUE )
    endif()
  endforeach()

  foreach(TPL_NAME ${${TRIBITS_PACKAGE}_TEST_REQUIRED_DEP_TPLS})
    if ( (NOT TPL_ENABLE_${TPL_NAME}) AND
      (NOT "${TPL_ENABLE_${TPL_NAME}}" STREQUAL "")
      )
      tribits_private_disable_tpl_required_package_enable(
        ${TPL_NAME}  ${TRIBITS_PACKAGE}  FALSE )
    endif()
  endforeach()

endmacro()


# Macro that disables all of the subpackages of a parent package.
#
macro(tribits_disable_parents_subpackages PARENT_PACKAGE_NAME)

  #message("TRIBITS_DISABLE_PARENTS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  #print_var(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

  if(NOT ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}
    AND NOT ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME} STREQUAL ""
    )

    set(SUBPACKAGE_IDX 0)
    foreach(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      set(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      set(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      #print_var(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
      if (NOT ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME} STREQUAL "OFF")
        set(ENABLE_BEING_DISABLED_VAR_NAME ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        message("-- "
          "Setting subpackage enable ${ENABLE_BEING_DISABLED_VAR_NAME}=OFF"
          " because parent package ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}=OFF")
        set(${ENABLE_BEING_DISABLED_VAR_NAME} OFF)
      endif()

      math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    endforeach()

  endif()

endmacro()


# Macro that enables all of the subpackages of a parent package.
#
macro(tribits_enable_parents_subpackages PARENT_PACKAGE_NAME)

  #message("TRIBITS_ENABLE_PARENTS_SUBPACKAGES: ${PARENT_PACKAGE_NAME}")

  #print_var(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

  if(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME})

    set(SUBPACKAGE_IDX 0)
    foreach(TRIBITS_SUBPACKAGE ${${PARENT_PACKAGE_NAME}_SUBPACKAGES})

      set(SUBPACKAGE_NAME ${TRIBITS_SUBPACKAGE})
      set(SUBPACKAGE_FULLNAME ${PARENT_PACKAGE_NAME}${TRIBITS_SUBPACKAGE})

      #print_var(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})

      if (NOT ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME} AND
        NOT "${${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}}" STREQUAL ""
        )
        # The subpackage is already disabled and is not just empty!
      elseif (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
        # The subpackage is already enabled so there is no reason to enable it!
      else()
        # The subpackage is not hard off or on so turn it on by default
        tribits_implicit_package_enable_is_allowed( "" ${SUBPACKAGE_FULLNAME}
          SUBPACKAGE_ALLOW_IMPLICIT_ENABLE)
        if (SUBPACKAGE_ALLOW_IMPLICIT_ENABLE)
          set(ENABLE_BEING_ENABLED_VAR_NAME ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
          message("-- "
            "Setting subpackage enable ${ENABLE_BEING_ENABLED_VAR_NAME}=ON"
            " because parent package ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}=ON")
          set(${ENABLE_BEING_ENABLED_VAR_NAME} ON)
        endif()
      endif()

      math(EXPR SUBPACKAGE_IDX "${SUBPACKAGE_IDX}+1")

    endforeach()

  endif()

endmacro()


# Macro that disables all forward packages that depend on the given packages
#
macro(tribits_disable_forward_required_dep_packages PACKAGE_NAME)

  #message("TRIBITS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES: ${PACKAGE_NAME}")

  if (
     (NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
     AND
     (NOT "${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}" STREQUAL "")
     )

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      tribits_private_disable_required_package_enables(${FWD_DEP_PKG} ${PACKAGE_NAME} TRUE)
    endforeach()

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      tribits_private_disable_optional_package_enables(${FWD_DEP_PKG} ${PACKAGE_NAME})
    endforeach()

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES})
      tribits_private_disable_required_package_enables(${FWD_DEP_PKG} ${PACKAGE_NAME} FALSE)
    endforeach()

  endif()

endmacro()


# Macro that prints out dependencies for a package
#
# Does not modify the global state.
#
macro(tribits_print_package_dependencies PACKAGE_NAME)

  set(PRINTED_VAR "")

  print_nonempty_var_with_spaces(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES PRINTED_VAR)

  if (${PROJECT_NAME}_DUMP_FORWARD_PACKAGE_DEPENDENCIES)
    print_nonempty_var_with_spaces(${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES
      PRINTED_VAR)
    print_nonempty_var_with_spaces(${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES
      PRINTED_VAR)
    print_nonempty_var_with_spaces(${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES
      PRINTED_VAR)
    print_nonempty_var_with_spaces(${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES
      PRINTED_VAR)
  endif()

  print_nonempty_var_with_spaces(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS PRINTED_VAR)
  print_nonempty_var_with_spaces(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS PRINTED_VAR)

  if (NOT PRINTED_VAR)
    message("-- ${PACKAGE_NAME}: No dependencies!")
  endif()

endmacro()


#
# Private helper macros
#


macro(tribits_private_add_optional_package_enable PACKAGE_NAME  OPTIONAL_DEP_PACKAGE
  TYPE  SET_AS_CACHE_IN
  )

  #message("\nPACKAGE_ARCH_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE}")

  if (SET_AS_CACHE_IN)

    multiline_set(DOCSTR
      "Enable optional ${TYPE} support in the package ${PACKAGE_NAME}"
      " for the package ${OPTIONAL_DEP_PACKAGE}."
      "  Set to 'ON', 'OFF', or leave empty"
      " to allow for other logic to decide."
      )

    set_cache_on_off_empty( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} ""
      ${DOCSTR} )

  else()

    if (NOT DEFINED ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
      set( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} "" )
    endif()

  endif()

endmacro()


macro(tribits_private_add_optional_tpl_enable PACKAGE_NAME OPTIONAL_DEP_TPL
  TYPE  SET_AS_CACHE_IN )

  if (SET_AS_CACHE_IN)

    multiline_set(DOCSTR
      "Enable optional ${TYPE} support in the package ${PACKAGE_NAME}"
      " for the TPL ${OPTIONAL_DEP_TPL}."
      "  Set to 'ON', 'OFF', or leave empty"
      " to allow for other logic to decide."
      )

    set_cache_on_off_empty( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} ""
      ${DOCSTR} )

  else()

    if (NOT DEFINED ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL})
      set( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} "" )
    endif()

  endif()

endmacro()


# Macro that sets cache vars for optional package interdependencies
#
# This also will set ${PACKAGE_NAME}_ENABLE_TESTS and
# ${PACKAGE_NAME}_ENABLE_EXAMPLES to empty non-cache vars
#
macro(tribits_set_up_optional_package_enables_and_cache_vars PACKAGE_NAME)

  #message("\nPACKAGE_ARCH_ADD_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  set(SET_AS_CACHE ${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}})

  if (SET_AS_CACHE)

    multiline_set(DOCSTR
      "Build tests for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set_cache_on_off_empty( ${PACKAGE_NAME}_ENABLE_TESTS "" ${DOCSTR} )

    multiline_set(DOCSTR
      "Build examples for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set_cache_on_off_empty( ${PACKAGE_NAME}_ENABLE_EXAMPLES "" ${DOCSTR} )

    multiline_set(DOCSTR
      "Build examples for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set( ${PACKAGE_NAME}_SKIP_CTEST_ADD_TEST
      "${${PROJECT_NAME}_SKIP_CTEST_ADD_TEST}" CACHE BOOL ${DOCSTR} )

  else()

    if (NOT DEFINED ${PACKAGE_NAME}_ENABLE_TESTS)
      set( ${PACKAGE_NAME}_ENABLE_TESTS "" )
    endif()
    if (NOT DEFINED ${PACKAGE_NAME}_ENABLE_EXAMPLES)
      set( ${PACKAGE_NAME}_ENABLE_EXAMPLES "" )
    endif()

  endif()

  foreach(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    tribits_private_add_optional_package_enable(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} "library" "${SET_AS_CACHE}" )
  endforeach()

  foreach(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
    tribits_private_add_optional_package_enable(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} "test" "${SET_AS_CACHE}" )
  endforeach()

  foreach(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
    tribits_private_add_optional_tpl_enable(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} "library" "${SET_AS_CACHE}" )
  endforeach()

  foreach(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
    tribits_private_add_optional_tpl_enable(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} "test" "${SET_AS_CACHE}" )
  endforeach()

endmacro()


# Macro that sets up the flat list of direct package dependencies and enabled
# package dependencies and sets ${packageName}_ENABLE_${depPkg} for LIB
# dependencies
#
# This makes it easy to just loop over all of the direct upstream dependencies
# for a package or just the enabled dependencies.
#
# NOTES:
#
#  * ${packageName}_LIB_DEFINED_DEPENDENCIES will be set regardless if
#    ${packageName} is enabled or not.
#
#  * ${packageName}_LIB_ENABLED_DEPENDENCIES is only set if ${packageName} is
#    enabled and will only contain the names of direct library upstream
#    internal and external packages ${depPkg} that are required or are
#    optional and ${packageName}_ENABLE_${depPkg} is set to ON.
#
#  * ${packageName}_TEST_DEFINED_DEPENDENCIES will be set regardless if
#    ${packageName} is enabled or not.
#
#  * ${packageName}_TEST_ENABLED_DEPENDENCIES is only set if ${packageName} is
#    enabled and will only contain the names of direct test/example upstream
#    internal and external packages ${depPkg} that are required or are
#    optional and ${packageName}_ENABLE_${depPkg} is set to ON.
#
#  * Sets ${packageName}_ENABLE_${depPkg}=ON for every required dep package
#    for LIB dependencies (but not TEST dependencies).  This allows looping
#    over just ${packageName}_LIB_DEFINED_DEPENDENCIES looking at
#    ${packageName}_ENABLE_${depPkg} to see if the package is enable or not.
#    This also includes special logic for required subpackages for parent
#    packages where only the shell of the parent package is enabled and not
#    all of its required subpackages are enabled.
#
macro(tribits_setup_direct_package_dependencies_lists_and_lib_required_enable_vars
    packageName
  )

  # LIB dependencies

  set(${packageName}_LIB_DEFINED_DEPENDENCIES "")
  set(${packageName}_LIB_ENABLED_DEPENDENCIES "")

  foreach(depPkg ${${packageName}_LIB_REQUIRED_DEP_PACKAGES})
    list(APPEND ${packageName}_LIB_DEFINED_DEPENDENCIES ${depPkg})
    if (${PROJECT_NAME}_ENABLE_${packageName} AND ${PROJECT_NAME}_ENABLE_${depPkg})
      set(${packageName}_ENABLE_${depPkg} ON)
      list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()
  # See below NOTE about required subpackage dependencies not being enabled in
  # some cases!

  foreach(depPkg ${${packageName}_LIB_OPTIONAL_DEP_PACKAGES})
    list(APPEND ${packageName}_LIB_DEFINED_DEPENDENCIES ${depPkg})
    if (${PROJECT_NAME}_ENABLE_${packageName} AND ${packageName}_ENABLE_${depPkg})
      list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  foreach(depPkg ${${packageName}_LIB_REQUIRED_DEP_TPLS})
    list(APPEND ${packageName}_LIB_DEFINED_DEPENDENCIES ${depPkg})
    if (${PROJECT_NAME}_ENABLE_${packageName})
      set(${packageName}_ENABLE_${depPkg} ON)
      list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  foreach(depPkg ${${packageName}_LIB_OPTIONAL_DEP_TPLS})
    list(APPEND ${packageName}_LIB_DEFINED_DEPENDENCIES ${depPkg})
    if (${PROJECT_NAME}_ENABLE_${packageName} AND ${packageName}_ENABLE_${depPkg})
      list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  # TEST dependencies

  set(${packageName}_TEST_DEFINED_DEPENDENCIES "")
  set(${packageName}_TEST_ENABLED_DEPENDENCIES "")

  if (${PROJECT_NAME}_ENABLE_${packageName}
      AND
      (${packageName}_ENABLE_TESTS OR ${packageName}_ENABLE_EXAMPLES)
    )
    set(enablePkgAndTestsOrExamples ON)
  else()
    set(enablePkgAndTestsOrExamples OFF)
  endif()

  foreach(depPkg ${${packageName}_TEST_REQUIRED_DEP_PACKAGES})
    list(APPEND ${packageName}_TEST_DEFINED_DEPENDENCIES ${depPkg})
    if (enablePkgAndTestsOrExamples)
      list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  foreach(depPkg ${${packageName}_TEST_OPTIONAL_DEP_PACKAGES})
    list(APPEND ${packageName}_TEST_DEFINED_DEPENDENCIES ${depPkg})
    if (enablePkgAndTestsOrExamples AND ${packageName}_ENABLE_${depPkg})
      list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  foreach(depPkg ${${packageName}_TEST_REQUIRED_DEP_TPLS})
    list(APPEND ${packageName}_TEST_DEFINED_DEPENDENCIES ${depPkg})
    if (enablePkgAndTestsOrExamples)
      list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

  foreach(depPkg ${${packageName}_TEST_OPTIONAL_DEP_TPLS})
    list(APPEND ${packageName}_TEST_DEFINED_DEPENDENCIES ${depPkg})
    if (enablePkgAndTestsOrExamples AND ${packageName}_ENABLE_${depPkg})
      list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
    endif()
  endforeach()

endmacro()
# NOTE: Above, a required dependency of an enabled package may not actually be
# enabled if it is a required subpackage of a parent package and the parent
# package was not actually enabled due to a dependency but the shell of the
# parent package was only enabled at the very end.  This is one of the more
# confusing aspects of the TriBITS dependency system.


# Function to print the direct package dependency lists
#
function(tribits_print_direct_package_dependencies_lists  packageName)
  set(PRINTED_VAR "")
  message("")
  print_nonempty_var_with_spaces(${packageName}_LIB_ENABLED_DEPENDENCIES PRINTED_VAR)
  print_var_with_spaces(${packageName}_LIB_DEFINED_DEPENDENCIES PRINTED_VAR)
  print_nonempty_var_with_spaces(${packageName}_TEST_ENABLED_DEPENDENCIES PRINTED_VAR)
  print_nonempty_var_with_spaces(${packageName}_TEST_DEFINED_DEPENDENCIES PRINTED_VAR)
endfunction()


#
# Private helper macros
#


# Enable optional intra-package support for enabled target package
# ${PACKAGE_NAME} (i.e. ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} is assumed to
# be TRUE before calling this macro.
#
macro(tribits_private_postprocess_optional_package_enable PACKAGE_NAME OPTIONAL_DEP_PACKAGE)

  #message("TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE}")
  #print_var(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
  #print_var(${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})

  if (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} AND ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
    message("-- " "NOTE:"
      " ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}"
      " is already set!")
  elseif ("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}" STREQUAL "")
    if (${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
      message("-- " "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON"
       " since ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON AND"
       " ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON")
      set(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} ON)
    else()
      message("-- " "NOT setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON"
       " since ${OPTIONAL_DEP_PACKAGE} is NOT enabled at this point!")
    endif()
  elseif (NOT "${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}" STREQUAL ""
    AND NOT ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}
    AND ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}
    )
    message("-- " "NOTE: ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}"
     " is already set so not enabling even though ${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=${${PROJECT_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}} is set!")
  endif()

  string(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  string(TOUPPER ${OPTIONAL_DEP_PACKAGE} OPTIONAL_DEP_PACKAGE_UPPER)
  set(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_PACKAGE_UPPER})

  if(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
    set(${MACRO_DEFINE_NAME} ON)
  else()
    set(${MACRO_DEFINE_NAME} OFF)
  endif()

endmacro()


# Enable optional intra-package support for enabled target package
# ${PACKAGE_NAME} (i.e. ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} is assumed to
# be TRUE before calling this macro.
#
macro(tribits_private_postprocess_optional_tpl_enable  PACKAGE_NAME  OPTIONAL_DEP_TPL)

  #message("TRIBITS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL}")

  if (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} AND TPL_ENABLE_${OPTIONAL_DEP_TPL})
    message("-- " "NOTE:"
      " ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}"
      " is already set!")
  elseif (
    (NOT TPL_ENABLE_${OPTIONAL_DEP_TPL})
    AND
    (NOT "${TPL_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL "")
    AND
    ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}
    )
    message(
      "\n***"
      "\n*** NOTE: Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=OFF"
      " which was ON since TPL_ENABLE_${OPTIONAL_DEP_TPL}=OFF"
      "\n***\n"
      )
    set(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} OFF)
  elseif ("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL ""
    AND TPL_ENABLE_${OPTIONAL_DEP_TPL}
    )
    message("-- " "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=ON"
      " since TPL_ENABLE_${OPTIONAL_DEP_TPL}=ON")
    set(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} ON)
  elseif (NOT "${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL ""
    AND NOT ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}
    AND TPL_ENABLE_${OPTIONAL_DEP_TPL}
    )
    message("-- " "NOTE: ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}"
      " is already set so not enabling even though TPL_ENABLE_${OPTIONAL_DEP_TPL}=${TPL_ENABLE_${OPTIONAL_DEP_TPL}} is set!")
  endif()

  string(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  string(TOUPPER ${OPTIONAL_DEP_TPL} OPTIONAL_DEP_TPL_UPPER)
  set(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_TPL_UPPER})

  if (${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL})
    set(${MACRO_DEFINE_NAME} ON)
  else()
    set(${MACRO_DEFINE_NAME} OFF)
  endif()

endmacro()


# Macro that post-processes optional dependencies after all other
# dependencies have been worked out
#
macro(tribits_postprocess_optional_package_enables PACKAGE_NAME)

  #message("\nPACKAGE_ARCH_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
      tribits_private_postprocess_optional_package_enable(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
    endforeach()

    foreach(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
      tribits_private_postprocess_optional_package_enable(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
    endforeach()

  endif()

endmacro()


# Macro that post-processes final package enables for packages with subpackage
# enables.
#
macro(tribits_postprocess_package_with_subpackages_enables  PACKAGE_NAME)
  #message("TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_ENABLES  '${PACKAGE_NAME}'")
  foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #print_var(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
    #print_var(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    #print_var(${SUBPACKAGE_FULLNAME}_ENABLE_TESTS)
    #print_var(${PACKAGE_NAME}_ENABLE_TESTS)
    #print_var(${SUBPACKAGE_FULLNAME}_ENABLE_EXAMPLES)
    #print_var(${PACKAGE_NAME}_ENABLE_EXAMPLES)
    if (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}
        AND NOT ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
      )
      message("-- "
        "Setting ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON"
        " because ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}=ON")
      set(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} ON)
      tribits_postprocess_package_with_subpackages_optional_subpackage_enables(
        ${PACKAGE_NAME})
      tribits_postprocess_package_with_subpackages_test_example_enables(
        ${PACKAGE_NAME}  TESTS)
      tribits_postprocess_package_with_subpackages_test_example_enables(
        ${PACKAGE_NAME}  EXAMPLES)
      # NOTE: We need to enable the parent package even if it was disabled by
      # some means before this because a subpackage is enabled.  But other
      # logic should ensure that the parent package is never disabled and a
      # subpackage is allowed to be enabled.
    endif()
  endforeach()
endmacro()


# Set <ParentPackage>_ENABLE_<SubPackage>=ON if not already enabled for all
# subpackages of a parent package.
#
macro(tribits_postprocess_package_with_subpackages_optional_subpackage_enables
    PACKAGE_NAME
  )
  #message("TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_TEST_EXAMPLE_ENABLES  '${PACKAGE_NAME}'")
  foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    if (${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}
      AND "${${PACKAGE_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}}" STREQUAL ""
      )
      message("-- "
        "Setting ${PACKAGE_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}=ON"
        " because ${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME}=ON")
      set(${PACKAGE_NAME}_ENABLE_${SUBPACKAGE_FULLNAME} ON)
    endif()
  endforeach()
endmacro()


# Set the parent package tests/examples enables if one subpackage is enabled
# and has its tests/examples
#
macro(tribits_postprocess_package_with_subpackages_test_example_enables
    PACKAGE_NAME  TESTS_OR_EXAMPLES
  )
  #message("TRIBITS_POSTPROCESS_PACKAGE_WITH_SUBPACKAGES_TEST_EXAMPLE_ENABLES  '${PACKAGE_NAME}'")
  foreach(TRIBITS_SUBPACKAGE ${${PACKAGE_NAME}_SUBPACKAGES})
    set(SUBPACKAGE_FULLNAME ${PACKAGE_NAME}${TRIBITS_SUBPACKAGE})
    #print_var(${PROJECT_NAME}_ENABLE_${SUBPACKAGE_FULLNAME})
    #print_var(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    #print_var(${SUBPACKAGE_FULLNAME}_ENABLE_${TESTS_OR_EXAMPLES})
    #print_var(${PACKAGE_NAME}_ENABLE_${TESTS_OR_EXAMPLES})
    #print_var(${SUBPACKAGE_FULLNAME}_ENABLE_EXAMPLES)
    #print_var(${PACKAGE_NAME}_ENABLE_EXAMPLES)
    if (${SUBPACKAGE_FULLNAME}_ENABLE_${TESTS_OR_EXAMPLES}
      AND "${${PACKAGE_NAME}_ENABLE_${TESTS_OR_EXAMPLES}}" STREQUAL ""
      )
      message("-- "
        "Setting ${PACKAGE_NAME}_ENABLE_${TESTS_OR_EXAMPLES}=ON"
        " because ${SUBPACKAGE_FULLNAME}_ENABLE_${TESTS_OR_EXAMPLES}=ON")
      set(${PACKAGE_NAME}_ENABLE_${TESTS_OR_EXAMPLES} ON)
    endif()
  endforeach()
endmacro()


# Post-processes optional package TPL based on if the TPL
# has been enabled or not
#
macro(tribits_postprocess_optional_tpl_enables PACKAGE_NAME)

  #message("\nPACKAGE_ARCH_ADD_OPTIONAL_TPL_ENABLES: ${PACKAGE_NAME}")
  
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
      tribits_private_postprocess_optional_tpl_enable(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
    endforeach()

    foreach(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
      tribits_private_postprocess_optional_tpl_enable(
        ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
    endforeach()

  endif()

endmacro()


# Set an individual package variable enable based on the global value
#
macro(tribits_set_all_packages_package_enable_variable   PACKAGE_ARCH_VAR   PACKAGE_VAR)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("")
    message("TRIBITS_SET_ALL_PACKAGES_PACKAGE_ENABLE_VARIABLE:")
    message("-- " "${PACKAGE_ARCH_VAR} = ${${PACKAGE_ARCH_VAR}}")
    message("-- " "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  endif()

  if ("${${PACKAGE_VAR}}" STREQUAL "")
    if (${PACKAGE_ARCH_VAR})
      message("-- " "Setting ${PACKAGE_VAR}=ON")
      set(${PACKAGE_VAR} ON)
    elseif (
      (NOT ${PACKAGE_ARCH_VAR})
      AND
      (NOT "${PACKAGE_ARCH_VAR}" STREQUAL "")
      )
      message("-- " "Setting ${PACKAGE_VAR}=OFF")
      set(${PACKAGE_VAR} OFF)
    else()
      if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        message("-- " "ELSE")
        # Otherwise, we will leave it up the the individual package
        # to decide?
      endif()
    endif()
  else()
    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message("-- " "${PACKAGE_VAR} NOT DEFAULT")
    endif()
  endif()

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("-- " "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  endif()

endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} based on
# ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
#
macro(tribits_apply_all_package_enables  PACKAGE_NAME)
  tribits_is_primary_meta_project_package(${PACKAGE_NAME}  PACKAGE_IS_PMPP)
  tribits_implicit_package_enable_is_allowed( "" ${PACKAGE_NAME}
    PROCESS_PACKAGE_ENABLE )
  if (PACKAGE_IS_PMPP  AND  PROCESS_PACKAGE_ENABLE)
    tribits_set_all_packages_package_enable_variable(
      ${PROJECT_NAME}_ENABLE_ALL_PACKAGES ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME} )
  endif()
endmacro()


# Macro used to set ${TRIBITS_PACKAGE)_ENABLE_TESTS and ${TRIBITS_PACKAGE)_ENABLE_EXAMPLES
# based on ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
#
macro(tribits_apply_test_example_enables PACKAGE_NAME)
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
    tribits_is_primary_meta_project_package(${PACKAGE_NAME}  PACKAGE_IS_PMPP)
    if (PACKAGE_IS_PMPP)
      tribits_set_all_packages_package_enable_variable(
        ${PROJECT_NAME}_ENABLE_TESTS ${PACKAGE_NAME}_ENABLE_TESTS )
      tribits_set_all_packages_package_enable_variable(
        ${PROJECT_NAME}_ENABLE_EXAMPLES ${PACKAGE_NAME}_ENABLE_EXAMPLES )
    endif()
  endif()
endmacro()


# Macro to disable ${PARENT_PACKAGE_NAME)_ENABLE_ENABLES by default if
# ${PARENT_PACKAGE_NAME)_ENABLE_TESTS is explicitly disabled.
#
macro(tribits_apply_package_examples_disable  PARENT_PACKAGE_NAME)
  if (NOT "${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}" STREQUAL ""
    AND NOT ${PARENT_PACKAGE_NAME}_ENABLE_TESTS
    AND "${${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
    )
    message("-- " "Setting"
      " ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES"
      "=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}"
      " because"
      " ${PARENT_PACKAGE_NAME}_ENABLE_TESTS"
      "=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}" )
     set(${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES ${${PARENT_PACKAGE_NAME}_ENABLE_TESTS})
  endif()
endmacro()
# NOTE: Above, the top-level package ${PARENT_PACKAGE_NAME} may not even be
# enabled yet when this gets called but its subpackages might and we need to
# process this default disable in case their are any enabled subpackages.


# Macro to disable ${TRIBITS_SUBPACKAGE)_ENABLE_TESTS and
# ${TRIBITS_SUBPACKAGE)_ENABLE_EXAMPLES based on
# ${TRIBITS_PARENTPACKAGE)_ENABLE_TESTS or
# ${TRIBITS_PARENTPACKAGE)_ENABLE_EXAMPLES
#
macro(tribits_apply_subpackage_tests_or_examples_disables  PARENT_PACKAGE_NAME
    TESTS_OR_EXAMPLES
  )
  set(parentPkgEnableVar ${PARENT_PACKAGE_NAME}_ENABLE_${TESTS_OR_EXAMPLES})
  if (NOT "${${parentPkgEnableVar}}" STREQUAL "" AND NOT ${parentPkgEnableVar})
    foreach(spkg IN LISTS ${PARENT_PACKAGE_NAME}_SUBPACKAGES)
      set(fullSpkgName ${PARENT_PACKAGE_NAME}${spkg})
      if (${PROJECT_NAME}_ENABLE_${fullSpkgName} AND NOT ${parentPkgEnableVar})
        if ("${${fullSpkgName}_ENABLE_${TESTS_OR_EXAMPLES}}" STREQUAL "")
          message("-- " "Setting"
            " ${fullSpkgName}_ENABLE_${TESTS_OR_EXAMPLES}=${${parentPkgEnableVar}}"
            " because parent package"
            " ${parentPkgEnableVar}=${${parentPkgEnableVar}}")
          set(${fullSpkgName}_ENABLE_${TESTS_OR_EXAMPLES} ${${parentPkgEnableVar}})
        endif()
      endif()
    endforeach()
  endif()
endmacro()


# Macro to enable ${TRIBITS_SUBPACKAGE)_ENABLE_TESTS and
# ${TRIBITS_SUBPACKAGE)_ENABLE_EXAMPLES based on
# ${TRIBITS_PARENTPACKAGE)_ENABLE_TESTS or
# ${TRIBITS_PARENTPACKAGE)_ENABLE_EXAMPLES
#
macro(tribits_apply_subpackage_tests_examples_enables  PARENT_PACKAGE_NAME)
  if ("${${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES}" STREQUAL ""
    AND ${PARENT_PACKAGE_NAME}_ENABLE_TESTS
    )
    message("-- " "Setting"
      " ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}"
      " because"
      " ${PARENT_PACKAGE_NAME}_ENABLE_TESTS=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}")
    set(${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES ${${PARENT_PACKAGE_NAME}_ENABLE_TESTS})
  endif()
  foreach(spkg IN LISTS ${PARENT_PACKAGE_NAME}_SUBPACKAGES)
    set(fullSpkgName ${PARENT_PACKAGE_NAME}${spkg})
    if (${PROJECT_NAME}_ENABLE_${fullSpkgName})
      if (${PARENT_PACKAGE_NAME}_ENABLE_TESTS)
        if ("${${fullSpkgName}_ENABLE_TESTS}" STREQUAL "")
          message("-- " "Setting"
            " ${fullSpkgName}_ENABLE_TESTS=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}"
            " because parent package"
            " ${PARENT_PACKAGE_NAME}_ENABLE_TESTS"
            "=${${PARENT_PACKAGE_NAME}_ENABLE_TESTS}")
          set(${fullSpkgName}_ENABLE_TESTS ${${PARENT_PACKAGE_NAME}_ENABLE_TESTS})
        endif()
      endif()
      if (${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES)
        if ("${${fullSpkgName}_ENABLE_EXAMPLES}" STREQUAL "")
          message("-- " "Setting"
            " ${fullSpkgName}_ENABLE_EXAMPLES=${${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES}"
            " because parent package"
            " ${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES"
            "=${${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES}")
          set(${fullSpkgName}_ENABLE_EXAMPLES ${${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES})
        endif()
      endif()
    endif()
  endforeach()
endmacro()
# NOTE: Above, the parent package may not actually be enabled yet
# (i.e. ${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME} my not be TRUE) if only
# subpackages needed to be enabled in the forward sweep but we want the tests
# and examples for subpackage to be enabled if
# ${PARENT_PACKAGE_NAME}_ENABLE_TESTS=ON or just examples i
# f${PARENT_PACKAGE_NAME}_ENABLE_EXAMPLES=ON


macro(tribits_private_enable_forward_package  FORWARD_DEP_PACKAGE_NAME  PACKAGE_NAME)
  tribits_implicit_package_enable_is_allowed( "" ${FORWARD_DEP_PACKAGE_NAME}
    ALLOW_PACKAGE_ENABLE )
  assert_defined(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
  if("${${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}}" STREQUAL ""
    AND ALLOW_PACKAGE_ENABLE
    )
    message("-- " "Setting ${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME}=ON"
      " because ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON")
    assert_defined(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
    set(${PROJECT_NAME}_ENABLE_${FORWARD_DEP_PACKAGE_NAME} ON)
  endif()
endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME)=ON for all
#
macro(tribits_enable_forward_lib_package_enables PACKAGE_NAME)

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      tribits_private_enable_forward_package(${FWD_DEP_PKG} ${PACKAGE_NAME})
    endforeach()

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      tribits_private_enable_forward_package(${FWD_DEP_PKG} ${PACKAGE_NAME})
    endforeach()

  endif()

endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME)=ON for all
# optional and required forward test dependencies of the package
# ${PACKAGE_NAME}
#
macro(tribits_enable_forward_test_package_enables PACKAGE_NAME)

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES})
      tribits_private_enable_forward_package(${FWD_DEP_PKG} ${PACKAGE_NAME})
    endforeach()

    foreach(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES})
      tribits_private_enable_forward_package(${FWD_DEP_PKG} ${PACKAGE_NAME})
    endforeach()

  endif()

endmacro()


#
# Private helper macros
#


macro(tribits_private_enable_dep_package  PACKAGE_NAME  DEP_PACKAGE_NAME
  OPTREQ_IN
  )

  #message("TRIBITS_PRIVATE_ENABLE_DEP_PACKAGE:  '${PACKAGE_NAME}'  '${DEP_PACKAGE_NAME}'  '${OPTREQ_IN}'")

  assert_defined(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})
  #print_var(${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})

    #message("The package is already enabled so there is nothing to enable!")

  elseif (${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME} STREQUAL "")

    set(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE "")

    if ("${OPTREQ_IN}" STREQUAL "REQUIRED")

      #message("Always enable the upstream dependency if it is required")

      message("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
        " because ${PACKAGE_NAME} has a required dependence on ${DEP_PACKAGE_NAME}")

      set(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)

    elseif (${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME})

      # Enable the upstream package if the user directly specified the
      # optional package enable regardless if it is PT or ST or even EX.

      message("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
        " because ${PACKAGE_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON")

      set(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)

    elseif (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)

      # Enable the package if there is an optional dependence and we are asked
      # to enabled optional dependencies.

      tribits_implicit_package_enable_is_allowed(${PACKAGE_NAME} ${DEP_PACKAGE_NAME}
        ALLOW_IMPLICIT_ENABLE)
      if (ALLOW_IMPLICIT_ENABLE)
        message("-- " "Setting ${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME}=ON"
          " because ${PACKAGE_NAME} has an optional dependence on ${DEP_PACKAGE_NAME}")
        set(TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE ON)
      endif()

    endif()

    # Enable the upstream package
    if (TRIBITS_PRIVATE_ENABLE_DEP_PACKAGES_ENABLE_PACKAGE)
      assert_defined(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME})
      set(${PROJECT_NAME}_ENABLE_${DEP_PACKAGE_NAME} ON)
    endif()

  endif()

endmacro()


macro(tribits_private_enable_dep_tpl  PACKAGE_NAME  DEP_TPL_NAME)
  assert_defined(TPL_ENABLE_${DEP_TPL_NAME})
  if(TPL_ENABLE_${DEP_TPL_NAME} STREQUAL "")
    message("-- " "Setting TPL_ENABLE_${DEP_TPL_NAME}=ON because"
      " it is required by the enabled package ${PACKAGE_NAME}")
    assert_defined(TPL_ENABLE_${DEP_TPL_NAME})
    set(TPL_ENABLE_${DEP_TPL_NAME} ON)
    set(TPL_${DEP_TPL_NAME}_ENABLING_PKG  ${PACKAGE_NAME})
  endif()
endmacro()


macro(tribits_private_enable_optional_dep_tpl PACKAGE_NAME DEP_TPL_NAME)
  #assert_defined(${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME})
  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
    AND ${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME}
    AND TPL_ENABLE_${DEP_TPL_NAME} STREQUAL ""
    )
    message("-- " "Setting TPL_ENABLE_${DEP_TPL_NAME}=ON because"
      " ${PACKAGE_NAME}_ENABLE_${DEP_TPL_NAME}=ON")
    assert_defined(TPL_ENABLE_${DEP_TPL_NAME})
    set(TPL_ENABLE_${DEP_TPL_NAME} ON)
  endif()
endmacro()


# Macro that enables the optional TPLs for given package
#
macro(tribits_enable_optional_tpls PACKAGE_NAME)

  #message("TRIBITS_ENABLE_OPTIONAL_TPLS: ${PACKAGE_NAME}")
  #message("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
      tribits_private_enable_optional_dep_tpl(${PACKAGE_NAME} ${DEP_TPL})
    endforeach()

    foreach(DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
      tribits_private_enable_optional_dep_tpl(${PACKAGE_NAME} ${DEP_TPL})
    endforeach()

  endif()

endmacro()


# Macro that enables upstream (required and optional) packages given package
#
# Here I have to enable the required packages too or the logic just does no
# work as expected.
#
macro(tribits_enable_upstream_packages  PACKAGE_NAME)

  #message("TRIBITS_ENABLE_UPSTREAM_PACKAGES: ${PACKAGE_NAME}")
  #message("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(DEP_PKG ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES})
      tribits_private_enable_dep_package(${PACKAGE_NAME} ${DEP_PKG} REQUIRED)
    endforeach()

    foreach(DEP_PKG ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
      tribits_private_enable_dep_package(${PACKAGE_NAME} ${DEP_PKG} OPTIONAL)
    endforeach()

    foreach(DEP_PKG ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES})
      tribits_private_enable_dep_package(${PACKAGE_NAME} ${DEP_PKG} REQUIRED)
    endforeach()

    foreach(DEP_PKG ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
      tribits_private_enable_dep_package(${PACKAGE_NAME} ${DEP_PKG} OPTIONAL)
    endforeach()

  endif()

endmacro()


# Macro that sets the required TPLs for given package
#
macro(tribits_enable_required_tpls PACKAGE_NAME)

  #message("PACKAGE_ARCH_ENABLE_REQUIRED_TPL_ENABLES: ${PACKAGE_NAME}")
  #message("-- " "${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=${${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}}")

  assert_defined(${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

  if (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})

    foreach(DEP_TPL ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS})
      tribits_private_enable_dep_tpl(${PACKAGE_NAME} ${DEP_TPL})
    endforeach()

    foreach(DEP_TPL ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS})
      tribits_private_enable_dep_tpl(${PACKAGE_NAME} ${DEP_TPL})
    endforeach()

  endif()

endmacro()


# @MACRO: tribits_adjust_package_enables()
#
# Usage:
#
#   tribits_adjust_package_enables()
#
# Macro that adjusts all of the package enables from what the user input to
# the final set that will be used to enable packages.
#
macro(tribits_adjust_package_enables)

  if (${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES)
    message("")
    message("Setting to empty '' all enabled packages on request ...")
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
      if (${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
        set_cache_on_off_empty(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} ""
          "Forced to empty '' by ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES=ON" FORCE)
        set(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE} "")
      endif()
      #print_var(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      # NOTE: Above, we don't want to set to empty those packages that have hard
      # disables because this will mess up the logic in later invocations.
    endforeach()
    advanced_set(${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
      "Forced to FALSE after use" FORCE)
  endif()

  #
  # A) Sweep forward through and apply all disables first!
  #

  tribits_get_nondisabled_list( ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_NOTDISABLED_PACKAGES "")

  message("")
  message("Disabling all packages that have a required dependency"
    " on disabled TPLs and optional package TPL support based on TPL_ENABLE_<TPL>=OFF ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_PACKAGES})
    tribits_disable_package_if_tpl_disabled(${TRIBITS_PACKAGE})
  endforeach()

  message("")
  message("Disabling subpackages for hard disables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=OFF ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
    tribits_disable_parents_subpackages(${TRIBITS_PACKAGE})
  endforeach()

  message("")
  message("Disabling forward required packages and optional intra-package"
    " support that have a dependency on disabled packages"
    " ${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
    tribits_disable_forward_required_dep_packages(${TRIBITS_PACKAGE})
  endforeach()

  tribits_get_nondisabled_list( ${PROJECT_NAME}_NOTDISABLED_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_NOTDISABLED_PACKAGES "")

  set(${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES
    "${${PROJECT_NAME}_NOTDISABLED_PACKAGES}")
  list(REVERSE ${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES)

  #
  # B) Apply all forward enables
  #

  message("")
  message("Enabling subpackages for hard enables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=ON ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_PACKAGES})
    tribits_enable_parents_subpackages(${TRIBITS_PACKAGE})
  endforeach()

  if (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    message("")
    message("Enabling all packages that are not currently disabled because of"
      " ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON"
      " (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})"
      " ...")
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_PACKAGES})
      tribits_apply_all_package_enables(${TRIBITS_PACKAGE})
    endforeach()
  endif()

  if (${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES)
    message("")
    message("Sweep forward enabling all forward library dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...")
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_NOTDISABLED_PACKAGES})
      tribits_enable_forward_lib_package_enables(${TRIBITS_PACKAGE})
    endforeach()
    message("")
    message("Sweep backward enabling all forward test dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...")
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES})
      tribits_enable_forward_test_package_enables(${TRIBITS_PACKAGE})
    endforeach()
    # NOTE: Above, we want to sweep backward to enable test-dependent packages
    # because we don't want to enable package Z just because package Y was enabled
    # because it had a test-only dependency on package X.  Sweeping backwards through
    # the packages makes sure this does not happen.
    set(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  endif()

  tribits_get_enabled_list( ${PROJECT_NAME}_NOTDISABLED_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES  "")

  #
  # C) Disable and enable tests for currently enabled packages
  #

  message("")
  message("Disabling subpackage tests/examples based on parent package tests/examples disables ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
    tribits_apply_package_examples_disable(${TRIBITS_PACKAGE} TESTS)
    tribits_apply_subpackage_tests_or_examples_disables(${TRIBITS_PACKAGE} TESTS)
    tribits_apply_subpackage_tests_or_examples_disables(${TRIBITS_PACKAGE} EXAMPLES)
  endforeach()

  if (${PROJECT_NAME}_ENABLE_TESTS OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    message("")
    message("Enabling all tests and/or examples that have not been"
      " explicitly disabled because ${PROJECT_NAME}_ENABLE_[TESTS,EXAMPLES]=ON ...")
    message("")
    foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
      tribits_apply_test_example_enables(${TRIBITS_PACKAGE})
    endforeach()
  endif()
  # NOTE: Above, we enable tests and examples here, before the remaining required
  # packages so that we don't enable tests that don't need to be enabled based
  # on the use of the option ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES.

  message("")
  message("Enabling subpackage tests/examples based on parent package tests/examples enables ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
    tribits_apply_subpackage_tests_examples_enables(${TRIBITS_PACKAGE})
  endforeach()
  # NOTE: We want to apply this logic here instead of later after the backward
  # sweep of package enables because we don't want to enable the
  # tests/examples for a subpackage if it is not needed based on set of
  # requested subpackages and packages to be enabled and the optional forward
  # sweep of downstream packages.

  #
  # D) Sweep backwards and enable upstream required and optional packages
  #

  if (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)
    set(EXTRA_MSG_STR " (and optional since ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON)")
  else()
    set(EXTRA_MSG_STR "")
  endif()

  message("")
  message("Enabling all required${EXTRA_MSG_STR} upstream packages for current set of"
    " enabled packages"
    " (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE=${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE})"
    " ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES})
    tribits_enable_upstream_packages(${TRIBITS_PACKAGE})
  endforeach()
  # NOTE: Above, we have to loop through the packages backward to enable all
  # the packages that feed into these packages.  This has to include *all*
  # upstream package enables including required packages, optional packages
  # (when ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES), and packages

  tribits_get_enabled_list( ${PROJECT_NAME}_NOTDISABLED_PACKAGES  ${PROJECT_NAME}
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES  "")

  message("")
  message("Enabling all optional intra-package enables <TRIBITS_PACKAGE>_ENABLE_<DEPPACKAGE>"
    " that are not currently disabled if both sets of packages are enabled ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    tribits_postprocess_optional_package_enables(${TRIBITS_PACKAGE})
  endforeach()

  #
  # E) Enable TPLs
  #

  message("")
  message("Enabling all remaining required TPLs for current set of"
    " enabled packages ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    tribits_enable_required_tpls(${TRIBITS_PACKAGE})
  endforeach()

  message("")
  message("Enabling all optional package TPL support"
    " <TRIBITS_PACKAGE>_ENABLE_<DEPTPL> not currently disabled for"
    " enabled TPLs ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    tribits_postprocess_optional_tpl_enables(${TRIBITS_PACKAGE})
  endforeach()

  message("")
  message("Enabling TPLs based on <TRIBITS_PACKAGE>_ENABLE_<TPL>=ON if TPL is not explicitly disabled ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    tribits_enable_optional_tpls(${TRIBITS_PACKAGE})
  endforeach()
  # NOTE: We need to do this after the above optional package TPL support
  # logic so that the TPL will be turned on for this package only as requested
  # in bug 4298.

  #
  # F) Set user cache variables for current set of enabled packages
  #

  message("")
  message("Set cache entries for optional packages/TPLs and tests/examples for packages actually enabled ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    tribits_set_up_optional_package_enables_and_cache_vars(${TRIBITS_PACKAGE})
  endforeach()

  #
  # G) Turn on parent packages where at least one subpackage has been enabled
  #

  message("")
  message("Enabling the shell of non-enabled parent packages (mostly for show) that have at least one subpackage enabled ...")
  message("")
  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES})
    tribits_postprocess_package_with_subpackages_enables(${TRIBITS_PACKAGE})
  endforeach()
  # NOTE: The above ensures that loops involving the parent package will
  # process the parent package but doing this last ensures that no downstream
  # dependencies will be enabled.

  tribits_set_up_enabled_lists_and_pkg_idx()

  #
  # H) Set up flat list of direct external and inner package dependencies (even
  # for non-enabled packages) and enabled package dependencies for enabled
  # packages
  #

  foreach(externalPkgName ${${PROJECT_NAME}_DEFINED_TPLS})
    tribits_extpkg_setup_enabled_dependencies(${externalPkgName})
    # ToDo: Assert that all of the listed dependencies in
    # ${externalPkgName}_LIB_ENABLED_DEPENDENCIES exist and are upstream from
    # ${externalPkgName}
  endforeach()

  foreach(TRIBITS_PACKAGE ${${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES})
    tribits_setup_direct_package_dependencies_lists_and_lib_required_enable_vars(
      ${TRIBITS_PACKAGE})
  endforeach()

  if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    message("\nDumping direct dependencies for each package ...")
    foreach(tribitsPkg  IN  LISTS  ${PROJECT_NAME}_DEFINED_TPLS
        ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
      )
      tribits_print_direct_package_dependencies_lists(${tribitsPkg})
    endforeach()
  endif()

endmacro()


# Function that sets up the full package dependencies for each enabled package
# including all of its indirect upstream package dependencies.
#
# This is needed in several different parts of the TriBITS implementation.
#
# ToDo: #63: Remove this function since we should not need a full list of
# direct and indirect package dependencies!
#
function(tribits_package_set_full_enabled_dep_packages  PACKAGE_NAME)

  set(PACKAGE_FULL_DEPS_LIST "")

  foreach(DEP_PKG ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES})
    if (${PROJECT_NAME}_ENABLE_${DEP_PKG})
      list(APPEND  PACKAGE_FULL_DEPS_LIST  ${DEP_PKG})
    endif()
    # NOTE: This if() should not be needed but this is a safeguard
  endforeach()

  foreach(DEP_PKG ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    if (${PACKAGE_NAME}_ENABLE_${DEP_PKG})
      list(APPEND  PACKAGE_FULL_DEPS_LIST  ${DEP_PKG})
    endif()
  endforeach()

  if(PACKAGE_FULL_DEPS_LIST)
    list(REMOVE_DUPLICATES  PACKAGE_FULL_DEPS_LIST)

    foreach(DEP_PACKAGE  ${PACKAGE_FULL_DEPS_LIST})
      list(APPEND PACKAGE_FULL_DEPS_LIST  ${${DEP_PACKAGE}_FULL_ENABLED_DEP_PACKAGES})
    endforeach()

    list(REMOVE_DUPLICATES PACKAGE_FULL_DEPS_LIST)
  endif()

  set(ORDERED_PACKAGE_FULL_DEPS_LIST "")

  foreach(DEP_PACKAGE  ${PACKAGE_FULL_DEPS_LIST})

    #print_var(${DEP_PACKAGE}_PKG_IDX)
    set(DEP_PACKAGE_VALUE  ${${DEP_PACKAGE}_PKG_IDX})

    set(SORTED_INDEX 0)
    set(INSERTED_DEP_PACKAGE FALSE)

    foreach(SORTED_PACKAGE  ${ORDERED_PACKAGE_FULL_DEPS_LIST})

      #print_var(${SORTED_PACKAGE}_PKG_IDX)
      set(SORTED_PACKAGE_VALUE  ${${SORTED_PACKAGE}_PKG_IDX})

      if (${DEP_PACKAGE_VALUE} GREATER ${SORTED_PACKAGE_VALUE})
        list(INSERT  ORDERED_PACKAGE_FULL_DEPS_LIST  ${SORTED_INDEX}  ${DEP_PACKAGE})
        set(INSERTED_DEP_PACKAGE TRUE)
        break()
      endif()

      math(EXPR SORTED_INDEX ${SORTED_INDEX}+1)

    endforeach()

    if(NOT INSERTED_DEP_PACKAGE)
      list(APPEND  ORDERED_PACKAGE_FULL_DEPS_LIST  ${DEP_PACKAGE})
    endif()

  endforeach()

  global_set(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES
    ${ORDERED_PACKAGE_FULL_DEPS_LIST})

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    print_var(${PACKAGE_NAME}_FULL_ENABLED_DEP_PACKAGES)
  endif()

endfunction()


# Function that creates enable-only dependency data-structures
#
# ToDo: #63: Remove this function since we should not need a full list of
# direct and indirect package dependencies!
#
function(tribits_set_up_enabled_only_dependencies)

  set(GENERATE_EXPORT_DEPENDENCIES ${${PROJECT_NAME}_GENERATE_EXPORT_FILE_DEPENDENCIES})
  set(lastExportTribitsPackage)

  if ("${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES}" STREQUAL ""
      AND NOT
      "${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}" STREQUAL ""
    )
    message(DEPRECATION
      "WARNING! The cache var"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES"
      "='${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}'"
      " is deprecated!  Please instead set"
      " ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES"
      "='${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES}'")
    set(${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES
      ${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_SE_PACKAGES} )
  endif()

  if (GENERATE_EXPORT_DEPENDENCIES
      AND ${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES
    )
    # Find the last enabled package for which an export file is requested.
    set(LAST_PKG_IDX -1)
    set(LAST_PKG)
    foreach(tribitsPkg ${${PROJECT_NAME}_GENERATE_EXPORT_FILES_FOR_ONLY_LISTED_PACKAGES})
      #print_var(tribitsPkg)
      set(PKG_IDX ${${tribitsPkg}_PKG_IDX})
      #print_var(PKG_IDX)
      if (PKG_IDX)
        # The listed package is enabled so we will consider it
        if (PKG_IDX GREATER ${LAST_PKG_IDX})
          set(LAST_PKG_IDX ${PKG_IDX})
          set(LAST_PKG ${tribitsPkg})
         #print_var(LAST_PKG_IDX)
         #print_var(LAST_PKG)
        endif()
      endif()
    endforeach()
    if (LAST_PKG)
      # At least one listed package was enabled
      set(lastExportTribitsPackage ${LAST_PKG})
    else()
      # None of the listed packages were enabled so don't bother generating
      # any export dependencies
      set(GENERATE_EXPORT_DEPENDENCIES FALSE)
    endif()

  endif()

  if (GENERATE_EXPORT_DEPENDENCIES)

    if (lastExportTribitsPackage)
      message("\nSetting up export dependencies up through ${lastExportTribitsPackage} ...\n")
    else()
      message("\nSetting up export dependencies for all enabled packages ...\n")
    endif()

    foreach(tribitsPackage ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
      tribits_package_set_full_enabled_dep_packages(${tribitsPackage})
      if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
        set(PRINTED_VAR FALSE)
        print_nonempty_var_with_spaces(${tribitsPackage}_FULL_ENABLED_DEP_PACKAGES
          PRINTED_VAR)
        if (NOT PRINTED_VAR)
          message("-- ${tribitsPackage}: No library dependencies!")
        endif()
      endif()
      if ("${lastExportTribitsPackage}" STREQUAL ${tribitsPackage})
        break()
      endif()
    endforeach()

  endif()

endfunction()
