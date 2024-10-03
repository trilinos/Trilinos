# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


include(TribitsProcessPackagesAndDirsLists)
include(TribitsAddOptionAndDefine)
include(TribitsGeneralMacros)
include(TribitsPrintEnabledPackagesLists)
include(TribitsPrintDependencyInfo)
include(TribitsPackageDependencies)
include(TribitsGetPackageEnableStatus)

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
include(TribitsCreateReverseList)


# NOTE: A nice way to view and navigate the contents of this file is to use
# the emacs 'occur' string:
#
#   "\(^##########\|^# .*-level\|^function\|^macro\)"


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
  tribits_unenable_enabled_packages()
  tribits_sweep_forward_apply_disables()
  tribits_sweep_forward_apply_enables()
  tribits_disable_and_enable_tests_and_examples()
  tribits_sweep_backward_enable_upstream_packages()
  tribits_set_cache_vars_for_current_enabled_packages()
  tribits_do_final_parent_packages_enables_for_subpackage_enables()
  tribits_adjust_internal_external_packages()
  tribits_setup_enabled_lists_and_pkg_idxs()
  tribits_setup_direct_packages_dependencies_lists_and_lib_required_enable_vars()
  tribits_print_direct_packages_dependencies_lists()
endmacro()



################################################################################
#
# First-level macros called directly from ``tribits_adjust_package_enables()``
#
# NOTE: In the below macros, local variables are prefixed by 'tap1_' in all of
# the macros() which will not clash because they are at the same level in the
# call stack (and are initialized in each macro).
# 
################################################################################


# @MACRO: tribits_unenable_enabled_packages()
#
# Macro to enable all unenabled packages if asked.
#
# See implementation for details.
#
macro(tribits_unenable_enabled_packages)
  if (${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES)
    message("\nSetting to empty '' all enabled packages because"
      " ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES="
      "'${${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES}'\n")
    foreach(tap1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES)
      if (${PROJECT_NAME}_ENABLE_${tap1_tribitsPkg})
        message("-- " "Setting ${PROJECT_NAME}_ENABLE_${tap1_tribitsPkg}=''")
        set_cache_on_off_empty(${PROJECT_NAME}_ENABLE_${tap1_tribitsPkg} ""
          "Forced to empty '' by ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES=ON" FORCE)
      endif()
      # NOTE: Above, we don't want to set to empty those packages that have hard
      # disables because this will mess up the logic in later invocations.
    endforeach()
    advanced_set(${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
      "Forced to FALSE after use" FORCE)
  endif()
endmacro()


# @MACRO: tribits_sweep_forward_apply_disables()
#
# Sweep forward and apply all disables in order first.
#
# See implementation for details.
#
macro(tribits_sweep_forward_apply_disables)

  message("\nDisabling all packages that have a required dependency"
    " on disabled TPLs and optional package TPL support based on TPL_ENABLE_<TPL>=OFF ...\n")

  message("\nDisabling subpackages for hard disables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=OFF ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES)
    tribits_disable_parents_subpackages(${tad1_tribitsPkg})
  endforeach()

  message("\nDisabling forward required packages and optional intra-package"
    " support that have a dependency on disabled packages"
    " ${PROJECT_NAME}_ENABLE_<TRIBITS_PACKAGE>=OFF (or"
    " TPL_ENABLE_<TRIBITS_EXTERNAL_PACKAGE>=OFF) ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_DEFINED_PACKAGES)
    tribits_disable_forward_required_dep_packages(${tad1_tribitsPkg})
  endforeach()

endmacro()


# @MACRO: tribits_sweep_forward_apply_enables()
#
# Updates the following variables in the local scope:
#
#   * ``${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES``
#   * `${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES`_
#
# See implementation for details.
#
macro(tribits_sweep_forward_apply_enables)

  tribits_get_sublist_nondisabled( ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES "")
  tribits_create_reverse_list(${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_REVERSE_NOTDISABLED_INTERNAL_PACKAGES)

  message("\nEnabling subpackages for hard enables of parent packages"
    " due to ${PROJECT_NAME}_ENABLE_<PARENT_PACKAGE>=ON ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES)
    tribits_enable_parents_subpackages(${tad1_tribitsPkg})
  endforeach()

  if (${PROJECT_NAME}_ENABLE_ALL_PACKAGES)
    message("\nEnabling all packages that are not currently disabled because of"
      " ${PROJECT_NAME}_ENABLE_ALL_PACKAGES=ON"
      " (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE="
      "${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}) ...\n")
    foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES)
      tribits_apply_all_package_enables(${tad1_tribitsPkg})
    endforeach()
  endif()

  if (${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES)
    message("\nSweep forward enabling all forward library dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...\n")
    foreach(tad1_tribitsPkg   IN LISTS  ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES)
      tribits_enable_forward_lib_package_enables(${tad1_tribitsPkg})
    endforeach()

    tribits_get_sublist_enabled( ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
      ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES "")
    tribits_create_reverse_list(${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES
      ${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES)

    message("\nSweep backward enabling all forward test dependent packages because"
      " ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON ...\n")
    foreach(tad1_tribitsPkg   IN LISTS
        ${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES
      )
      tribits_enable_forward_test_package_enables(${tad1_tribitsPkg})
    endforeach()
    # NOTE: Above, we want to sweep backward to enable test-dependent packages
    # because we don't want to enable package Z just because package Y was enabled
    # because it had a test-only dependency on package X.  Sweeping backwards through
    # the packages makes sure this does not happen.
    set(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  endif()

  tribits_get_sublist_enabled( ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES  "")

endmacro()
# NOTE: Above, we are sweeping over *all* of the not-disabled packages listed
# in ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES, including those package that
# might have <Package>_PACKAGE_BUILD_STATUS=EXTERNAL.  That makes sense
# because these are TriBITS (or TriBITS-compliant) packages so we should
# assume that all of their downstream packages, whether internal or external,
# should be enabled as well.  If we find this is not the desirable behavior,
# then we can change this later.


# @MACRO: tribits_disable_and_enable_tests_and_examples()
#
# Adjust test and example enables based on different criteria.
#
# See implementation for details.
#
macro(tribits_disable_and_enable_tests_and_examples)

  message("\nDisabling subpackage tests/examples based on parent package"
    " tests/examples disables ...\n")
  foreach(tad1_tribitsPkg  IN LISTS
      ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    )
    tribits_apply_package_examples_disable(${tad1_tribitsPkg} TESTS)
    tribits_apply_subpackage_tests_or_examples_disables(${tad1_tribitsPkg} TESTS)
    tribits_apply_subpackage_tests_or_examples_disables(${tad1_tribitsPkg} EXAMPLES)
  endforeach()

  if (${PROJECT_NAME}_ENABLE_TESTS  OR  ${PROJECT_NAME}_ENABLE_EXAMPLES)
    message("\nEnabling all tests and/or examples that have not been"
      " explicitly disabled because ${PROJECT_NAME}_ENABLE_[TESTS,EXAMPLES]=ON ...\n")
    foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES)
      tribits_apply_test_example_enables(${tad1_tribitsPkg})
    endforeach()
  endif()
  # NOTE: Above, we enable tests and examples here, before the remaining required
  # packages so that we don't enable tests that don't need to be enabled based
  # on the use of the option ${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES.

  message("\nEnabling subpackage tests/examples based on parent package tests/examples enables ...\n")
  foreach(tad1_tribitsPkg  IN LISTS
      ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    )
    tribits_apply_subpackage_tests_examples_enables(${tad1_tribitsPkg})
  endforeach()
  # NOTE: We want to apply this logic here instead of later after the backward
  # sweep of package enables because we don't want to enable the
  # tests/examples for a subpackage if it is not needed based on set of
  # requested subpackages and packages to be enabled and the optional forward
  # sweep of downstream packages.

endmacro()


# @MACRO: tribits_sweep_backward_enable_upstream_packages()
#
# Sweep backwards and enable required (and optional) upstream packages.
#
# This sets the final value for:
#
# * `${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES`_
#
# See implementation for details.
#
macro(tribits_sweep_backward_enable_upstream_packages)

  set(tap1_extraMsgStr "")
  if (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)
    set(tap1_extraMsgStr
       " (and optional since ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES=ON)")
  endif()

  tribits_get_sublist_nondisabled( ${PROJECT_NAME}_DEFINED_PACKAGES
    ${PROJECT_NAME}_NOTDISABLED_PACKAGES "")
  tribits_create_reverse_list(${PROJECT_NAME}_NOTDISABLED_PACKAGES
    ${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES)

  message("\nEnabling all required${tap1_extraMsgStr} upstream packages for current"
    " set of enabled packages (${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE="
    "${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}) ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_REVERSE_NOTDISABLED_PACKAGES)
    tribits_enable_upstream_packages(${tad1_tribitsPkg})
  endforeach()
  # NOTE: Above, we have to loop through the packages backward to enable all
  # the packages that this package depends on.  This has to include *all*
  # upstream package enables including required packages, and optional
  # packages (when ${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES).

  tribits_get_sublist_enabled( ${PROJECT_NAME}_NOTDISABLED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES  "")

  message("\nEnabling all optional intra-package enables"
    " <TRIBITS_PACKAGE>_ENABLE_<DEPPACKAGE> that are not currently disabled"
    " if both sets of packages are enabled ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES)
    tribits_postprocess_optional_package_enables(${tad1_tribitsPkg})
  endforeach()

endmacro()


# @MACRO: tribits_set_cache_vars_for_current_enabled_packages()
#
macro(tribits_set_cache_vars_for_current_enabled_packages)
  message("\nSet cache entries for optional packages/TPLs and tests/examples for packages actually enabled ...\n")
  foreach(tad1_tribitsPkg  IN LISTS  ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES)
    tribits_setup_optional_package_enables_and_cache_vars(${tad1_tribitsPkg})
  endforeach()
endmacro()


# @MACRO: tribits_do_final_parent_packages_enables_for_subpackage_enables()
#
macro(tribits_do_final_parent_packages_enables_for_subpackage_enables)

  message("\nEnabling the shell of non-enabled parent packages (mostly for show)"
    " that have at least one subpackage enabled ...\n")
  foreach(tad1_tribitsPkg  IN LISTS
      ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    )
    tribits_enable_parent_package_for_subpackage_enables(${tad1_tribitsPkg})
  endforeach()
  # NOTE: The above ensures that loops involving the parent package will
  # process the parent package.  But this has to be done after all of the
  # other enable/disable logic to ensure no downstream dependencies will be
  # enabled based on this.

endmacro()


# @MACRO: tribits_adjust_internal_external_packages()
#
# Macro to adjust the set of internal and external packages by changing
# `<Package>_PACKAGE_BUILD_STATUS`.
#
# This macro sweeps backwards over the dependency graph setting.
#
# NOTE: This is called **after** all of the logic that determines what
# packages are enabled or disabled.  This is because we don't want to change
# the enable/disable logic when one or more initially internal packages are
# made external.  We are going to just assume that if an initially internal
# package that is declared external should be disabled, then the user will
# need to make that decision explicitly.
#
macro(tribits_adjust_internal_external_packages)
  tribits_create_reverse_list(${PROJECT_NAME}_DEFINED_PACKAGES
    ${PROJECT_NAME}_REVERSE_DEFINED_PACKAGES)
  message("\nAdjust the set of internal and external packages:\n")
  foreach(packageName  IN LISTS  ${PROJECT_NAME}_REVERSE_DEFINED_PACKAGES)
    tribits_set_package_and_related_upstream_packages_to_external(${packageName})
  endforeach()
endmacro()


# Macro that sets up the basic lists of enabled packages and packages.
#
macro(tribits_setup_enabled_lists_and_pkg_idxs)

  # ${PROJECT_NAME}_ENABLED_PACKAGES
  tribits_get_sublist_enabled(
    ${PROJECT_NAME}_DEFINED_INTERNAL_TOPLEVEL_PACKAGES
    ${PROJECT_NAME}_ENABLED_INTERNAL_TOPLEVEL_PACKAGES
    ${PROJECT_NAME}_NUM_ENABLED_INTERNAL_TOPLEVEL_PACKAGES)

  # ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES
  tribits_get_sublist_enabled( ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES
    ${PROJECT_NAME}_NUM_ENABLED_INTERNAL_PACKAGES)

  # ${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES
  unset(${PROJECT_NAME}_REVERSE_ENABLED_INTERNAL_PACKAGES) # Wipe out temp value

  # Set ${tribitsPackage}_PKG_IDX for each enabled package
  set(pkgIdx 0)
  foreach(tribitsPackage ${${PROJECT_NAME}_ENABLED_INTERNAL_PACKAGES})
    set(${tribitsPackage}_PKG_IDX ${pkgIdx})
    math(EXPR  pkgIdx  "${pkgIdx} + 1")
  endforeach()

endmacro()


# @MACRO: tribits_setup_direct_packages_dependencies_lists_and_lib_required_enable_vars()
#
# Set up flat list of direct external and inner package dependencies (even for
# non-enabled packages) and enabled package dependencies for enabled packages.
#
macro(tribits_setup_direct_packages_dependencies_lists_and_lib_required_enable_vars)

  foreach(tad1_externalPkgName  IN LISTS  ${PROJECT_NAME}_DEFINED_TPLS)
    tribits_extpkg_setup_enabled_dependencies(${tad1_externalPkgName})
    # ToDo: Assert that all of the listed dependencies in
    # ${tad1_externalPkgName}_LIB_ENABLED_DEPENDENCIES exist and are
    # upstream from ${tad1_externalPkgName}
  endforeach()

  foreach(tad1_internalPkgName  IN LISTS  ${PROJECT_NAME}_DEFINED_INTERNAL_PACKAGES)
    tribits_setup_enabled_dependencies_lists_and_enable_vars(${tad1_internalPkgName})
  endforeach()

endmacro()


# @MACRO: tribits_print_direct_packages_dependencies_lists()
#
macro(tribits_print_direct_packages_dependencies_lists)
  if (${PROJECT_NAME}_DUMP_PACKAGE_DEPENDENCIES)
    message("\nDumping direct enabled dependencies for each package ...")
    foreach(tad1_tribitsPkg  IN  LISTS  ${PROJECT_NAME}_DEFINED_PACKAGES)
      tribits_print_direct_package_dependencies_lists(${tad1_tribitsPkg})
    endforeach()
  endif()
endmacro()



################################################################################
#
# Second-level macros called indirectly from ``tribits_adjust_package_enables()``
#
################################################################################


# Macro that disables all of the subpackages of a parent package.
#
macro(tribits_disable_parents_subpackages  parentPackageName)

  if(NOT ${PROJECT_NAME}_ENABLE_${parentPackageName}
      AND (NOT ${PROJECT_NAME}_ENABLE_${parentPackageName} STREQUAL "")
    )

    foreach(tap2_subPkgName   IN LISTS  ${parentPackageName}_SUBPACKAGES)

      set(subpkgFullName ${parentPackageName}${tap2_subPkgName})

      if (NOT ${PROJECT_NAME}_ENABLE_${subpkgFullName} STREQUAL "OFF")
        set(packageBeingDisabledVarName ${PROJECT_NAME}_ENABLE_${subpkgFullName})
        message("-- "
          "Setting subpackage enable ${packageBeingDisabledVarName}=OFF"
          " because parent package ${PROJECT_NAME}_ENABLE_${parentPackageName}=OFF")
        set(${packageBeingDisabledVarName} OFF)
      endif()

    endforeach()

  endif()

endmacro()


# Macro that disables forward package that depends on the passed-in package
#
macro(tribits_disable_forward_required_dep_packages  packageName)
  tribits_get_package_enable_status(${packageName}  packageEnable  "")
  if ((NOT packageEnable) AND (NOT "${packageEnable}" STREQUAL ""))
    foreach(fwdDepPkg  IN LISTS  ${packageName}_FORWARD_LIB_DEFINED_DEPENDENCIES)
      if (${fwdDepPkg}_LIB_DEP_REQUIRED_${packageName})
        tribits_private_disable_required_package_enables(${fwdDepPkg}
          ${packageName} TRUE)
      else()
        tribits_private_disable_optional_package_enables(${fwdDepPkg}
          ${packageName} TRUE)
      endif()
    endforeach()
    foreach(fwdDepPkg  IN LISTS  ${packageName}_FORWARD_TEST_DEFINED_DEPENDENCIES)
      if (${fwdDepPkg}_TEST_DEP_REQUIRED_${packageName})
        tribits_private_disable_required_package_enables(${fwdDepPkg}
          ${packageName} FALSE)
      endif()
    endforeach()
  endif()
endmacro()


# Macro that enables all of the subpackages of a parent package.
#
macro(tribits_enable_parents_subpackages  parentPackageName)

  if(${PROJECT_NAME}_ENABLE_${parentPackageName})

    set(subpkgIdx 0)
    foreach(tap2_subPkgName  IN LISTS  ${parentPackageName}_SUBPACKAGES)

      set(subpkgFullName ${parentPackageName}${tap2_subPkgName})

      if (NOT ${PROJECT_NAME}_ENABLE_${subpkgFullName} AND
        NOT "${${PROJECT_NAME}_ENABLE_${subpkgFullName}}" STREQUAL ""
        )
        # The subpackage is already disabled and is not just empty!
      elseif (${PROJECT_NAME}_ENABLE_${subpkgFullName})
        # The subpackage is already enabled so there is no reason to enable it!
      else()
        # The subpackage is not hard off or on so turn it on by default
        tribits_implicit_package_enable_is_allowed( "" ${subpkgFullName}
          subpkgAllowImplicitEnable)
        if (subpkgAllowImplicitEnable)
          set(enableVarName ${PROJECT_NAME}_ENABLE_${subpkgFullName})
          message("-- "
            "Setting subpackage enable ${enableVarName}=ON"
            " because parent package ${PROJECT_NAME}_ENABLE_${parentPackageName}="
	    "${${PROJECT_NAME}_ENABLE_${parentPackageName}}")
          set(${enableVarName} ON)
        endif()
      endif()

      math(EXPR  subpkgIdx "${subpkgIdx}+1")

    endforeach()

  endif()

endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_<packageName> based on
# ${PROJECT_NAME}_ENABLE_ALL_PACKAGES
#
macro(tribits_apply_all_package_enables  packageName)
  tribits_is_primary_meta_project_package(${packageName}  packageIsPmpp)
  tribits_implicit_package_enable_is_allowed( "" ${packageName}
    processThisPackageEnable )
  if (packageIsPmpp  AND  processThisPackageEnable)
    tribits_set_package_enable_based_on_project_enable(
      ${PROJECT_NAME}_ENABLE_ALL_PACKAGES  ${PROJECT_NAME}_ENABLE_${packageName} )
  endif()
endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME}=ON for all
# forward required or optional LIB dependency on ${packageName}
#
macro(tribits_enable_forward_lib_package_enables  packageName)
  if (${PROJECT_NAME}_ENABLE_${packageName})
    foreach(fwdDepPkgName  IN LISTS  ${packageName}_FORWARD_LIB_DEFINED_DEPENDENCIES)
      tribits_private_enable_forward_package(${fwdDepPkgName} ${packageName})
    endforeach()
  endif()
endmacro()


# Macro used to set ${PROJECT_NAME}_ENABLE_${FWD_PACKAGE_NAME}=ON for all
# forward required or optional TEST dependency on ${packageName}
#
macro(tribits_enable_forward_test_package_enables  packageName)
  foreach(fwdDepPkgName  IN LISTS  ${packageName}_FORWARD_TEST_DEFINED_DEPENDENCIES)
    tribits_private_enable_forward_package(${fwdDepPkgName} ${packageName})
  endforeach()
endmacro()
# NOTE: The above macro does not need to check if ${packageName} is enabled
# because it will only be called for packages that are enabled already.  This
# not only improves performance but it also aids in testing.


# Macro to disable ${parentPackageName)_ENABLE_ENABLES by default if
# ${parentPackageName)_ENABLE_TESTS is explicitly disabled.
#
macro(tribits_apply_package_examples_disable  parentPackageName)
  if ( (NOT ${parentPackageName}_ENABLE_TESTS)
      AND (NOT "${${parentPackageName}_ENABLE_TESTS}" STREQUAL "")
      AND ("${${parentPackageName}_ENABLE_EXAMPLES}" STREQUAL "")
    )
    message("-- " "Setting"
      " ${parentPackageName}_ENABLE_EXAMPLES"
      "=${${parentPackageName}_ENABLE_TESTS}"
      " because"
      " ${parentPackageName}_ENABLE_TESTS"
      "=${${parentPackageName}_ENABLE_TESTS}" )
     set(${parentPackageName}_ENABLE_EXAMPLES ${${parentPackageName}_ENABLE_TESTS})
  endif()
endmacro()
# NOTE: Above, the top-level package ${parentPackageName} may not even be
# enabled yet when this gets called but its subpackages might and we need to
# process this default disable in case their are any enabled subpackages.


# Macro to disable <tribitsSubPackage>_ENABLE_TESTS and
# <tribitsSubPackage>_ENABLE_EXAMPLES based on
# <parentPackageName>_ENABLE_TESTS or <parentPackageName>_ENABLE_EXAMPLES
#
macro(tribits_apply_subpackage_tests_or_examples_disables  parentPackageName
    testsOrExamples
  )
  set(parentPkgEnableVar ${parentPackageName}_ENABLE_${testsOrExamples})
  if ((NOT ${parentPkgEnableVar}) AND (NOT "${${parentPkgEnableVar}}" STREQUAL ""))
    foreach(subpkgName  IN LISTS  ${parentPackageName}_SUBPACKAGES)
      set(fullSpkgName ${parentPackageName}${subpkgName})
      if (${PROJECT_NAME}_ENABLE_${fullSpkgName})
        if ("${${fullSpkgName}_ENABLE_${testsOrExamples}}" STREQUAL "")
          message("-- " "Setting"
            " ${fullSpkgName}_ENABLE_${testsOrExamples}=${${parentPkgEnableVar}}"
            " because parent package"
            " ${parentPkgEnableVar}=${${parentPkgEnableVar}}")
          set(${fullSpkgName}_ENABLE_${testsOrExamples} ${${parentPkgEnableVar}})
        endif()
      endif()
    endforeach()
  endif()
endmacro()


# Macro used to set <packageName>_ENABLE_TESTS and
# <packageName>_ENABLE_EXAMPLES based on ${PROJECT_NAME}_ENABLE_TESTS and
# ${PROJECT_NAME}_ENABLE_EXAMPLES
#
macro(tribits_apply_test_example_enables  packageName)
  if (${PROJECT_NAME}_ENABLE_${packageName})
    tribits_is_primary_meta_project_package(${packageName}  packageIsPmmp)
    if (packageIsPmmp)
      tribits_set_package_enable_based_on_project_enable_on(
        ${PROJECT_NAME}_ENABLE_TESTS  ${packageName}_ENABLE_TESTS )
      tribits_set_package_enable_based_on_project_enable_on(
        ${PROJECT_NAME}_ENABLE_EXAMPLES  ${packageName}_ENABLE_EXAMPLES )
    endif()
  endif()
endmacro()


# Macro to enable subpackage tests and examples based on parent package tests
# and examples enables
#
macro(tribits_apply_subpackage_tests_examples_enables  parentPackageName)
  # Set default for ${parentPackageName}_ENABLE_EXAMPLES=OFF if tests disabled
  if ( ("${${parentPackageName}_ENABLE_EXAMPLES}" STREQUAL "")
      AND ${parentPackageName}_ENABLE_TESTS
    )
    message("-- " "Setting"
      " ${parentPackageName}_ENABLE_EXAMPLES=${${parentPackageName}_ENABLE_TESTS}"
      " because"
      " ${parentPackageName}_ENABLE_TESTS=${${parentPackageName}_ENABLE_TESTS}")
    set(${parentPackageName}_ENABLE_EXAMPLES ${${parentPackageName}_ENABLE_TESTS})
  endif()
  # Set defaults for <fullSubpackageName>_ENABLE_[TESTS|EXAMPLES]
  set(parentEnableExamples ${${parentPackageName}_ENABLE_EXAMPLES})
  set(parentEnableTests ${${parentPackageName}_ENABLE_TESTS})
  foreach(subpkgName  IN LISTS  ${parentPackageName}_SUBPACKAGES)
    set(fullSpkgName ${parentPackageName}${subpkgName})
    if (${PROJECT_NAME}_ENABLE_${fullSpkgName})
      if (parentEnableTests AND ("${${fullSpkgName}_ENABLE_TESTS}" STREQUAL ""))
        message("-- " "Setting"
          " ${fullSpkgName}_ENABLE_TESTS=${parentEnableTests}"
          " because parent package"
          " ${parentPackageName}_ENABLE_TESTS=${parentEnableTests}")
        set(${fullSpkgName}_ENABLE_TESTS ${parentEnableTests})
      endif()
      if (parentEnableExamples AND ("${${fullSpkgName}_ENABLE_EXAMPLES}" STREQUAL ""))
        message("-- " "Setting"
          " ${fullSpkgName}_ENABLE_EXAMPLES=${parentEnableExamples}"
          " because parent package"
          " ${parentPackageName}_ENABLE_EXAMPLES=${parentEnableExamples}")
        set(${fullSpkgName}_ENABLE_EXAMPLES ${parentEnableExamples})
      endif()
    endif()
  endforeach()
endmacro()
# NOTE: Above, the parent package may not actually be enabled yet
# (i.e. ${PROJECT_NAME}_ENABLE_${parentPackageName} my not be TRUE) if only
# subpackages needed to be enabled in the forward sweep but we want the tests
# and examples for a subpackage to be enabled if the tests and examples for
# the parent package, respectfully, are enabled.


# Macro that enables upstream (required and optional) packages for a given
# package
#
macro(tribits_enable_upstream_packages  packageName)

  tribits_get_package_enable_status(${packageName}  packageEnable  packageEnableVar)

  if (packageEnable)

    foreach(depPkg  IN LISTS  ${packageName}_LIB_DEFINED_DEPENDENCIES)
      tribits_private_enable_dep_package(${packageName}  ${depPkg}  LIB)
    endforeach()

    foreach(depPkg  IN LISTS  ${packageName}_TEST_DEFINED_DEPENDENCIES)
      tribits_private_enable_dep_package(${packageName}  ${depPkg}  TEST)
    endforeach()

  endif()

endmacro()
# NOTE: The above macro has a defect.  It is enabling upstream test dependent
# packages even if tests nor examples are not enabled (see
# TriBITSPub/TriBITS#56).  But fixing this will break backward compatibility
# and therefore require upgrading the packages that currently only work
# correctly because of this defect.


# Macro that post-processes optional dependencies after all other
# dependencies have been worked out
#
macro(tribits_postprocess_optional_package_enables  packageName)

  if (${PROJECT_NAME}_ENABLE_${packageName})

    foreach(depPkg ${${packageName}_LIB_DEFINED_DEPENDENCIES})
      tribits_private_postprocess_optional_package_enable(
        ${packageName} ${depPkg} )
    endforeach()

    foreach(depPkg ${${packageName}_TEST_DEFINED_DEPENDENCIES})
      tribits_private_postprocess_optional_package_enable(
        ${packageName} ${depPkg} )
    endforeach()

  endif()

endmacro()
# NOTE: Above, it is harmless to process required dependencies as well so we
# leave of the if () statement based on
# ${packageName}_[LIB|TEST]_DEP_REQUIRED_${depPkg}.


# Macro that sets cache vars for optional package interdependencies
#
# This also will set ${PACKAGE_NAME}_ENABLE_TESTS and
# ${PACKAGE_NAME}_ENABLE_EXAMPLES to empty non-cache vars
#
macro(tribits_setup_optional_package_enables_and_cache_vars  packageName)

  assert_defined(${PROJECT_NAME}_ENABLE_${packageName})
  set(setAsCacheVar ${${PROJECT_NAME}_ENABLE_${packageName}})

  if (setAsCacheVar)

    multiline_set(docStr
      "Build tests for the package ${packageName}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set_cache_on_off_empty( ${packageName}_ENABLE_TESTS "" ${docStr} )

    multiline_set(docStr
      "Build examples for the package ${packageName}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set_cache_on_off_empty( ${packageName}_ENABLE_EXAMPLES "" ${docStr} )

    multiline_set(docStr
      "Build examples for the package ${packageName}.  Set to 'ON', 'OFF', or leave empty ''"
       " to allow for other logic to decide."
       )
    set( ${packageName}_SKIP_CTEST_ADD_TEST
      "${${PROJECT_NAME}_SKIP_CTEST_ADD_TEST}" CACHE BOOL ${docStr} )

  else()

    if (NOT DEFINED ${packageName}_ENABLE_TESTS)
      set( ${packageName}_ENABLE_TESTS "" )
    endif()
    if (NOT DEFINED ${packageName}_ENABLE_EXAMPLES)
      set( ${packageName}_ENABLE_EXAMPLES "" )
    endif()

  endif()

  foreach(optDepPkg ${${packageName}_LIB_DEFINED_DEPENDENCIES})
    tribits_private_add_optional_package_enable(
      ${packageName} ${optDepPkg} "library" "${setAsCacheVar}" )
  endforeach()

  foreach(optDepPkg ${${packageName}_TEST_DEFINED_DEPENDENCIES})
    tribits_private_add_optional_package_enable(
      ${packageName} ${optDepPkg} "test" "${setAsCacheVar}" )
  endforeach()

endmacro()


# Macro that enables a the top-level parent package enable if any of its
# subpackages are enabled.
#
macro(tribits_enable_parent_package_for_subpackage_enables  toplevelPackageName)
  foreach(tap2_subPkgName  IN LISTS  ${toplevelPackageName}_SUBPACKAGES)
    set(subpkgFullName ${toplevelPackageName}${tap2_subPkgName})
    tribits_get_package_enable_status(${toplevelPackageName}
      toplevelPackageEnable  toplevelPackageEnableVarName)
    tribits_get_package_enable_status(${subpkgFullName}
      subpkgEnable  subpkgEnableVarName)
    if (subpkgEnable AND (NOT toplevelPackageEnable))
      message("-- "
        "Setting ${toplevelPackageEnableVarName}=ON because"
	" ${subpkgEnableVarName}=${subpkgEnable}")
      set(${toplevelPackageEnableVarName} ON)
      tribits_set_parent_package_subpackage_enable_for_enabled_subpackages(
        ${toplevelPackageName})
      tribits_set_parent_package_test_example_enable_for_enabled_subpackages(
        ${toplevelPackageName}  TESTS)
      tribits_set_parent_package_test_example_enable_for_enabled_subpackages(
        ${toplevelPackageName}  EXAMPLES)
      # NOTE: Above, we need to enable the parent package even if it was
      # disabled by some means before this.  (There are use cases where the
      # parent package my be disabled but that may not trigger the disable of
      # subpackages of that package.)  Other logic should ensure that the
      # parent package is never explicitly disabled and a subpackage is
      # allowed to be enabled.
      break() # We only need trigger above code for single enabled subpackage!
    endif()
  endforeach()
endmacro()


# Macro that sets $``{packageName}_PACKAGE_BUILD_STATUS=EXTERNAL and direct
# upstream dependent packages ``${depPkg}_PACKAGE_BUILD_STATUS=EXTERNAL`` if
# required.
#
# * Set an INTERNAL package as EXTERNAL if ``TPL_ENABLE_${packageName}`` is
#  ``TRUE``
#
# * Set a top-level package to EXTERNAL if any of its subpackages are EXTERNAL
#   (or ``TPL_ENABLE_${subpkgFullName}`` is TRUE).
#
# This macro must be called in a backward/reverse sweep over the dependency
# graph to work correctly.  It sets upstream packages as EXTERNAL if the given
# downstream package ``<packageName>`` is EXTERNAL.  There is also special
# logic for handling a package that has subpackages.  That is, if any
# subpackage of ``<packageName>`` is determined is EXTERNAL, then the parent
# package is set to EXTERNAL and all of the other subpackages in parent
# package ``<packageName>`` are set as EXTERNAL as well.  (We don't allow a
# subset of subpackages in a parent package to be EXTERNAL and the other
# subpackages to be INTERNAL.  That would be way too complicated to implement
# and be way too confusing for implementors and users.)
#
macro(tribits_set_package_and_related_upstream_packages_to_external  packageName)

  tribits_set_parent_package_external_if_subpackage_external(${packageName}
    subpackageTriggeredParentPackageExternal)

  tribits_set_package_to_external_if_requested_by_user(${packageName})

  tribits_set_upstream_dep_packages_as_external(${packageName}
    ${subpackageTriggeredParentPackageExternal})

  tribits_mark_package_as_upstream_of_tribits_compliant_external_package(${packageName})

endmacro()
# NOTE: In the above macro, if ${packageName} is made EXTERNAL because it one
# of its subpackages is considered EXTERNAL, then the loop over all of the
# package's upstream dependencies listed in
# ${packageName}_LIB_DEFINED_DEPENDENCIES results in all of the other
# subpackages in that package to also be treated as external.  (I.e., there is
# no need for a special loop over subpackages for the parent package.  And the
# only direct dependencies for a parent package should be its subpackages.  If
# that is not the case, then, technically, this would be setting more packages
# as EXTERNAL than need to be.)


# Macro that sets up the flat list of direct package dependencies and enabled
# package dependencies and sets ${packageName}_ENABLE_${depPkg} for LIB
# dependencies
#
# This makes it easy to just loop over all of the direct upstream dependencies
# for a package or just the enabled dependencies.
#
# NOTES:
#
#  * ${packageName}_LIB_ENABLED_DEPENDENCIES is only set if ${packageName} is
#    enabled and will only contain the names of direct library upstream
#    internal and external packages ${depPkg} that are required or are
#    optional and ${packageName}_ENABLE_${depPkg} is set to ON.
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
macro(tribits_setup_enabled_dependencies_lists_and_enable_vars  packageName)

  tribits_get_package_enable_status(${packageName} packageEnable packageEnableVar)

  set(${packageName}_LIB_ENABLED_DEPENDENCIES "")
  foreach(depPkg IN LISTS ${packageName}_LIB_DEFINED_DEPENDENCIES)
    tribits_get_package_enable_status(${depPkg} depPkgEnable depPkgEnableVar)
    if (${packageName}_LIB_DEP_REQUIRED_${depPkg})
      if (packageEnable AND depPkgEnable)
        set(${packageName}_ENABLE_${depPkg} ON)
        # See below NOTE about required subpackage dependencies not being
        # enabled in some cases
        list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
      endif()
    else()
      if (packageEnable AND ${packageName}_ENABLE_${depPkg})
        list(APPEND ${packageName}_LIB_ENABLED_DEPENDENCIES ${depPkg})
      endif()
    endif()
  endforeach()

  set(${packageName}_TEST_ENABLED_DEPENDENCIES "")
  if (packageEnable AND (${packageName}_ENABLE_TESTS OR ${packageName}_ENABLE_EXAMPLES))
    foreach(depPkg IN LISTS ${packageName}_TEST_DEFINED_DEPENDENCIES)
      if (${packageName}_TEST_DEP_REQUIRED_${depPkg})
        list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
      else()
        if (${packageName}_ENABLE_${depPkg})
          list(APPEND ${packageName}_TEST_ENABLED_DEPENDENCIES ${depPkg})
        endif()
      endif()
    endforeach()
  endif()

endmacro()
# NOTE: Above, a required dependency of an enabled package may not actually be
# enabled if the upstream depPkg is a required subpackage of a parent package
# and the parent package was not actually enabled due to a dependency, but
# instead, only the shell of the parent package was enabled at the very end
# (in tribits_do_final_parent_packages_enables_for_subpackage_enables()).
# This is one of the more confusing aspects of the TriBITS dependency system.


# Function to print the direct package dependency lists
#
function(tribits_print_direct_package_dependencies_lists  packageName)
  message("")
  set(printedVar "")
  tribits_print_nonempty_package_deps_list(${packageName}  LIB  ENABLED  printedVar)
  tribits_print_nonempty_package_deps_list(${packageName}  TEST  ENABLED  printedVar)
  if (NOT printedVar)
    message("-- ${packageName}: No enabled dependencies!")
  endif()
endfunction()



################################################################################
#
# Third and lower-level macros called indirectly from
# ``tribits_adjust_package_enables()``
#
################################################################################


# Only turn off ${fwdDepPkgName} libraries or test/examples if it is currently
# enabled or could be enabled
#
macro(tribits_private_disable_required_package_enables
    fwdDepPkgName  packageName  libraryDep
  )
  tribits_get_package_enable_status(${fwdDepPkgName}  ""  fwdDepPkgEnableVarName)
  if (${fwdDepPkgEnableVarName} OR "${${fwdDepPkgEnableVarName}}" STREQUAL "")
    if ("${libraryDep}" STREQUAL "TRUE")
      tribits_private_print_disable_required_package_enable(
        ${packageName}  ${fwdDepPkgEnableVarName}
        ${fwdDepPkgName} "library" )
      set(${fwdDepPkgEnableVarName}  OFF)
    else()
      set(depTypeStr "test/example")
      if (${fwdDepPkgName}_ENABLE_TESTS
        OR "${${fwdDepPkgName}_ENABLE_TESTS}" STREQUAL ""
        )
        tribits_private_print_disable_required_package_enable(
          ${packageName} ${fwdDepPkgName}_ENABLE_TESTS
          ${fwdDepPkgName} "${depTypeStr}" )
        set(${fwdDepPkgName}_ENABLE_TESTS OFF)
      endif()

      if (${fwdDepPkgName}_ENABLE_EXAMPLES
        OR "${${fwdDepPkgName}_ENABLE_EXAMPLES}" STREQUAL ""
        )
        tribits_private_print_disable_required_package_enable(
          ${packageName} ${fwdDepPkgName}_ENABLE_EXAMPLES
          ${fwdDepPkgName} "${depTypeStr}" )
        set(${fwdDepPkgName}_ENABLE_EXAMPLES OFF)
      endif()
    endif()
  endif()
endmacro()

function(tribits_private_print_disable_required_package_enable
    packageName  packageEnableSomethingVarName  fwdDepPkgName
    depTypeStr
  )
  tribits_private_print_disable(
    ${packageEnableSomethingVarName} ${fwdDepPkgName}
    "${depTypeStr}" "package" ${packageName} )
endfunction()


function(tribits_private_print_disable
  packageBeingDisabledVarName  packageWithSomethingBeingDisabledName
  depTypeStr  thingBeingDisabledType  thingBeingDisabledName
  )
  if (${packageBeingDisabledVarName})
    if (${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES)
      message(
        " ***\n"
        " *** NOTE: Setting ${packageBeingDisabledVarName}=OFF"
        " which was '${${packageBeingDisabledVarName}}' because"
        " ${packageWithSomethingBeingDisabledName} has"
        " a required ${depTypeStr} dependence on disabled"
        " ${thingBeingDisabledType} ${thingBeingDisabledName}"
        " but ${PROJECT_NAME}_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON!\n"
        " ***\n"
        )
    else()
      message(FATAL_ERROR
        " ***\n"
        " *** ERROR: Setting ${packageBeingDisabledVarName}=OFF"
        " which was '${${packageBeingDisabledVarName}}' because"
        " ${packageWithSomethingBeingDisabledName} has"
        " a required ${depTypeStr} dependence on disabled"
        " ${thingBeingDisabledType} ${thingBeingDisabledName}!\n"
        " ***\n"
        )
    endif()
  else()
    message("-- "
      "Setting ${packageBeingDisabledVarName}=OFF"
      " because ${packageWithSomethingBeingDisabledName} has a required ${depTypeStr}"
      " dependence on disabled ${thingBeingDisabledType} ${thingBeingDisabledName}")
  endif()
endfunction()


macro(tribits_private_disable_optional_package_enables  fwdDepPkgName  packageName)

  if (${fwdDepPkgName}_ENABLE_${packageName}
      OR "${${fwdDepPkgName}_ENABLE_${packageName}}" STREQUAL ""
    )
    # Always disable the conditional enable but only print the message if the
    # package is enabled or if a disable overrides an enable
    if (${PROJECT_NAME}_ENABLE_${fwdDepPkgName})
      if (${fwdDepPkgName}_ENABLE_${packageName})  # is explicitly enabled already!
        message("-- "
          "NOTE: Setting ${fwdDepPkgName}_ENABLE_${packageName}=OFF"
          " which was ${${fwdDepPkgName}_ENABLE_${packageName}}"
          " because ${fwdDepPkgName} has an optional library dependence"
          " on disabled package ${packageName}")
      else()  # Not explicitly set
        message("-- "
          "Setting ${fwdDepPkgName}_ENABLE_${packageName}=OFF"
          " because ${fwdDepPkgName} has an optional library dependence"
          " on disabled package ${packageName}")
      endif()
    endif()
    if (${fwdDepPkgName}_ENABLE_${packageName}
        AND (NOT ${PROJECT_NAME}_ENABLE_${packageName})
        AND (NOT "${${PROJECT_NAME}_ENABLE_${packageName}}" STREQUAL "")
      )
      message("-- " "NOTE: ${fwdDepPkgName}_ENABLE_${packageName}="
        "${${fwdDepPkgName}_ENABLE_${packageName}} but"
        " ${PROJECT_NAME}_ENABLE_${packageName}="
        "${${PROJECT_NAME}_ENABLE_${packageName}} is set.  Setting"
        " ${fwdDepPkgName}_ENABLE_${packageName}=OFF!")
    endif()
    set(${fwdDepPkgName}_ENABLE_${packageName} OFF)
  endif()

endmacro()


# Set an individual package variable enable variable (to ON or OFF) based on a
# global enable value
#
macro(tribits_set_package_enable_based_on_project_enable  projectEnableVar
    packageEnableVar
  )

  if ("${${packageEnableVar}}" STREQUAL "")
    if (${projectEnableVar})
      message("-- " "Setting ${packageEnableVar}=ON")
      set(${packageEnableVar} ON)
    elseif ( (NOT ${projectEnableVar})
        AND (NOT "${projectEnableVar}" STREQUAL "")
      )
      message("-- " "Setting ${packageEnableVar}=OFF")
      set(${packageEnableVar} OFF)
    else()
      # Otherwise, we will leave it up the the individual package
      # to decide?
    endif()
  else()
    # "${packageEnableVar} not at the default empty ''
  endif()

endmacro()


# Set an individual package test or examples enable to on only if global
# enable var is on
#
macro(tribits_set_package_enable_based_on_project_enable_on  projectEnableVar
    packageEnableVar
  )
  if (("${${packageEnableVar}}" STREQUAL "") AND ${projectEnableVar})
    message("-- " "Setting ${packageEnableVar}=ON")
    set(${packageEnableVar} ON)
  endif()
endmacro()


macro(tribits_private_enable_dep_package  packageName  depPkgName  libOrTest)
  tribits_get_package_enable_status(${depPkgName}  depPkgEnable  depPkgEnableVar)
  if (depPkgEnable)
    #message("The package is already enabled so there is nothing to enable!")
  elseif ("${depPkgEnable}"  STREQUAL  "")
    set(tpedp_enableDepPkg "")
    if (${packageName}_${libOrTest}_DEP_REQUIRED_${depPkgName}
        AND (NOT "${${packageName}_SOURCE_DIR}" STREQUAL "")
      )
      message("-- " "Setting ${depPkgEnableVar}=ON"
        " because ${packageName} has a required dependence on ${depPkgName}")
      set(tpedp_enableDepPkg  ON)
    elseif (${packageName}_ENABLE_${depPkgName})
      # Enable the upstream package if the user directly specified the
      # optional package enable regardless if it is PT or ST or even EX.
      message("-- " "Setting ${depPkgEnableVar}=ON"
        " because ${packageName}_ENABLE_${depPkgName}=ON")
      set(tpedp_enableDepPkg  ON)
    elseif (${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES)
      # Enable the package if there is an optional dependence and we are asked
      # to enabled optional dependencies.
      tribits_implicit_package_enable_is_allowed(${packageName}  ${depPkgName}
        allowImplicitEnable)
      if (allowImplicitEnable)
        message("-- " "Setting ${depPkgEnableVar}=ON"
          " because ${packageName} has an optional dependence on ${depPkgName}")
        set(tpedp_enableDepPkg  ON)
      endif()
    endif()
    # Enable the upstream package
    if (tpedp_enableDepPkg)
      set(${depPkgEnableVar}  ON)
      set(${depPkgName}_ENABLING_PKG  ${packageName})
    endif()
  endif()
endmacro()


# Enable optional intra-package support for enabled target package
# ${packageName} (i.e. ${PROJECT_NAME}_ENABLE_${packageName} is assumed to
# be TRUE before calling this macro.
#
macro(tribits_private_postprocess_optional_package_enable  packageName  optDepPkg)

  tribits_get_package_enable_status(${optDepPkg}  optDepPkgEnable  optDepPkgEnableVar)
  tribits_get_package_enable_status(${packageName} packageEnable  packageEnableVar)

  if (${packageName}_ENABLE_${optDepPkg}  AND  optDepPkgEnable)
    message("-- " "NOTE:"
      " ${packageName}_ENABLE_${optDepPkg}=${${packageName}_ENABLE_${optDepPkg}}"
      " is already set!")
  elseif ("${${packageName}_ENABLE_${optDepPkg}}" STREQUAL "")
    if (optDepPkgEnable)
      message("-- " "Setting ${packageName}_ENABLE_${optDepPkg}=ON"
       " since ${packageEnableVar}=ON AND"
       " ${optDepPkgEnableVar}=ON")
      set(${packageName}_ENABLE_${optDepPkg}  ON)
    else()
      message("-- " "NOT setting ${packageName}_ENABLE_${optDepPkg}=ON"
       " since ${optDepPkg} is NOT enabled at this point!")
    endif()
  elseif ((NOT "${${packageName}_ENABLE_${optDepPkg}}" STREQUAL "")
      AND (NOT ${packageName}_ENABLE_${optDepPkg})
      AND optDepPkgEnable
    )
    message("-- " "NOTE: ${packageName}_ENABLE_${optDepPkg}="
      "${${packageName}_ENABLE_${optDepPkg}} is already set so not enabling even"
      " though ${optDepPkgEnableVar}="
      "${optDepPkgEnable} is set!")
  endif()

  string(TOUPPER  ${packageName}  packageName_UPPER)
  string(TOUPPER  ${optDepPkg}  optDepPkg_UPPER)
  set(macroDefineName  HAVE_${packageName_UPPER}_${optDepPkg_UPPER})

  if(${packageName}_ENABLE_${optDepPkg})
    set(${macroDefineName} ON)
  else()
    set(${macroDefineName} OFF)
  endif()

endmacro()


macro(tribits_private_enable_forward_package  fwdDepPkgName  packageName)
  tribits_implicit_package_enable_is_allowed( "" ${fwdDepPkgName} allowFwdDepPkgEnable)
  if ("${${PROJECT_NAME}_ENABLE_${fwdDepPkgName}}" STREQUAL "" AND allowFwdDepPkgEnable)
    message("-- " "Setting ${PROJECT_NAME}_ENABLE_${fwdDepPkgName}=ON"
      " because ${PROJECT_NAME}_ENABLE_${packageName}=ON")
    set(${PROJECT_NAME}_ENABLE_${fwdDepPkgName}  ON)
  endif()
endmacro()


macro(tribits_private_add_optional_package_enable  packageName  optionalDepPkgName
  libraryOrTest  setAsCacheVar
  )

  if (setAsCacheVar)

    multiline_set(docStr
      "Enable optional ${libraryOrTest} support in the package ${packageName}"
      " for the package ${optionalDepPkgName}."
      "  Set to 'ON', 'OFF', or leave empty"
      " to allow for other logic to decide."
      )

    set_cache_on_off_empty( ${packageName}_ENABLE_${optionalDepPkgName} ""
      ${docStr} )

  else()

    if (NOT DEFINED ${packageName}_ENABLE_${optionalDepPkgName})
      set(${packageName}_ENABLE_${optionalDepPkgName} "" )
    endif()

  endif()

endmacro()


# Set the parent package as EXTERNAL if any of its subpackages are being
# treated as EXTERNAL and return bool if that is the case.
#
# On output, ``<parentPackageName>_PACKAGE_BUILD_STATUS`` will have been set
# to ``EXTERNAL`` if the returned var
# ``subpackageTriggeredParentPackageExternalOut`` is ``TRUE``.
#
macro(tribits_set_parent_package_external_if_subpackage_external  parentPackageName
    subpackageTriggeredParentPackageExternalOut
  )

  set(subpackageTriggeredParentPackageExternal FALSE)

  if (NOT ${parentPackageName}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")
    foreach (tap2_subPkgName  IN LISTS  ${parentPackageName}_SUBPACKAGES)
      set(subpkgFullName ${parentPackageName}${tap2_subPkgName})
      tribits_package_is_external(${subpkgFullName} subpkgIsExternal)
      if (subpkgIsExternal)
	set(becauseMsgStr "subpackage ${subpkgFullName} being treated as EXTERNAL")
        if (TPL_ENABLE_${subpkgFullName})
          string(APPEND becauseMsgStr
	    " (TPL_ENABLE_${subpkgFullName}=${TPL_ENABLE_${subpkgFullName}})")
        endif()
        tribits_set_internal_package_to_external(${parentPackageName} "${becauseMsgStr}")
        set(subpackageTriggeredParentPackageExternal TRUE)
        break()
      endif()
    endforeach()
  endif()

  set(${subpackageTriggeredParentPackageExternalOut}
    ${subpackageTriggeredParentPackageExternal})

endmacro()


# Macro that sets ``<packageName>_PACKAGE_BUILD_STATUS`` to ``EXTERNAL`` if
# requested by the user (e.g. by setting ``TPL_ENABLE_<packageName>``).
#
macro(tribits_set_package_to_external_if_requested_by_user  packageName)
  if (NOT ${packageName}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")
    if (TPL_ENABLE_${packageName})
      tribits_set_internal_package_to_external(${packageName}
        "TPL_ENABLE_${packageName}=${TPL_ENABLE_${packageName}}")
    endif()
  endif()
endmacro()


# Macro that sets all of the direct upstream dependent packages as EXTERNAL if
# they should be.
#
# NOTE: We only bother setting upstream dependencies as EXTERNAL if the
# package ``<packageName>`` is enabled or if it is a parent package (enabled
# or not) where at least one of its subpackages is enabled and is being
# treated as EXTERNAL.  (In the latter case, any subpackages that are not
# enabled will still be set as EXTERNAL.)
#
macro(tribits_set_upstream_dep_packages_as_external  packageName
    subpackageTriggeredParentPackageExternal
  )

  tribits_get_package_enable_status(${packageName} packageEnable "")

  if (${packageName}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")
    foreach(depPkg  IN LISTS  ${packageName}_LIB_DEFINED_DEPENDENCIES)
      if ((NOT ${depPkg}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL") AND
          (subpackageTriggeredParentPackageExternal OR packageEnable)
        )
        tribits_set_internal_package_to_external(${depPkg}
          "downstream package ${packageName} being treated as EXTERNAL")
      endif()
    endforeach()
  endif()

endmacro()


# Mark a package as being upstream of a TriBITS-compliant external package
# (and therefore should be processed by a downstream TriBITS-complaint
# package)
#
# NOTE: These packages are initially marked by setting
# ``${packageName}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE=TRUE``.
# If this packages are not actually defined by a downstream TriBITS-compliant
# external package, then this variable will be set to ``FALSE`` later.
#
macro(tribits_mark_package_as_upstream_of_tribits_compliant_external_package  packageName)

  set_default(${packageName}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE FALSE)

  tribits_get_package_enable_status(${packageName} packageEnable "")

  if (${packageName}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")

    foreach(fwdDepPkg  IN LISTS  ${packageName}_FORWARD_LIB_DEFINED_DEPENDENCIES)

      if((${fwdDepPkg}_IS_TRIBITS_COMPLIANT
            OR ${fwdDepPkg}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE)
          AND (${fwdDepPkg}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")
        )
        tribits_get_package_enable_status(${fwdDepPkg} fwdDepPkgEnable "")
        if (${fwdDepPkg}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE)
          set(directOrIndirectStr "indirectly")
          set(downstreamPkgStr "")
        else()
          set(directOrIndirectStr "directly")
          set(downstreamPkgStr " ${fwdDepPkg}")
        endif()
        if (packageEnable AND (NOT ${packageName}_IS_TRIBITS_COMPLIANT))
          message("-- "
            "NOTE: ${packageName} is ${directOrIndirectStr} upstream from a"
            " TriBITS-compliant external package${downstreamPkgStr}")
        endif()
        set(${packageName}_PROCESSED_BY_DOWNSTREAM_TRIBITS_EXTERNAL_PACKAGE  TRUE)
        break()
      endif()

    endforeach()

  endif()

endmacro()



# Macro that sets ``<ParentPackage>_ENABLE_<SubPackage>=ON`` if not already
# enabled for all enabled subpackages of a parent package.
#
macro(tribits_set_parent_package_subpackage_enable_for_enabled_subpackages
    toplevelPackageName
  )
  foreach(tap3_subPkg   IN LISTS  ${toplevelPackageName}_SUBPACKAGES)
    set(subpkgFullName ${toplevelPackageName}${tap3_subPkg})
    if (${PROJECT_NAME}_ENABLE_${subpkgFullName}
      AND "${${toplevelPackageName}_ENABLE_${subpkgFullName}}" STREQUAL ""
      )
      message("-- "
        "Setting ${toplevelPackageName}_ENABLE_${subpkgFullName}=ON"
        " because ${PROJECT_NAME}_ENABLE_${subpkgFullName}=ON")
      set(${toplevelPackageName}_ENABLE_${subpkgFullName} ON)
    endif()
  endforeach()
endmacro()


# Macro that sets ``<ParentPacakge>_ENABLE_[TESTS|EXAMPLES]=ON`` if subpackage
# is enabled and has its tests/examples are enabled.
#
macro(tribits_set_parent_package_test_example_enable_for_enabled_subpackages
    toplevelPackageName  testOrExamples
  )
  foreach(tap3_subPkg  IN LISTS  ${toplevelPackageName}_SUBPACKAGES)
    set(subpkgFullName ${toplevelPackageName}${tap3_subPkg})
    if (${subpkgFullName}_ENABLE_${testOrExamples}
      AND "${${toplevelPackageName}_ENABLE_${testOrExamples}}" STREQUAL ""
      )
      message("-- "
        "Setting ${toplevelPackageName}_ENABLE_${testOrExamples}=ON"
        " because ${subpkgFullName}_ENABLE_${testOrExamples}=ON")
      set(${toplevelPackageName}_ENABLE_${testOrExamples} ON)
    endif()
  endforeach()
endmacro()


# Macro that sets an internal package to EXTERNAL and print so and why (if the
# package is actually enabled)
#
# Usage::
#
#   tribits_set_internal_package_to_external(<depPkgName> "<becauseMsg1>"
#     "<becauseMsg2>" ...)
#
# This always sets ``<depPkgName>_PACKAGE_BUILD_STATUS=EXTERNAL`` but only
# prints the message if ``<depPkgName>`` is enabled.
#
macro(tribits_set_internal_package_to_external  depPkgName)
  if (NOT ${depPkgName}_INTERNAL_PACKAGE_ALREADY_SET_EXTERNAL)
    tribits_get_package_enable_status(${depPkgName} depPkgEnable "")
    if (depPkgEnable)
      message("-- "
         "Treating internal package ${depPkgName} as EXTERNAL because"
         " " ${ARGN})
    endif()
    set(${depPkgName}_PACKAGE_BUILD_STATUS  EXTERNAL)
    set(${depPkgName}_FINDMOD  TRIBITS_PKG)
    set(${depPkgName}_INTERNAL_PACKAGE_ALREADY_SET_EXTERNAL TRUE)
  endif()
endmacro()


# Function to return if <packageName> is to be treated as an EXTERNAL package
# in processing of the package
#
function(tribits_package_is_external  packageName  packageIsExternalOut)
  if (TPL_ENABLE_${packageName})
    set(packageIsExternal TRUE)
  elseif (${packageName}_PACKAGE_BUILD_STATUS STREQUAL "EXTERNAL")
    set(packageIsExternal TRUE)
  else()
    set(packageIsExternal FALSE)
  endif()
  set(${packageIsExternalOut} ${packageIsExternal} PARENT_SCOPE)
endfunction()


# LocalWords: tribits TriBITS foreach endmacro endfunction
