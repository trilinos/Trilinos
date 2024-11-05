# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsPackageMacros)
include(TribitsReportInvalidTribitsUsage)


# @MACRO: tribits_subpackage()
#
# Forward declare a `TriBITS Subpackage`_ called at the top of the
# subpackage's `<packageDir>/<spkgDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   tribits_subpackage(<spkgName>)
#
# Once called, the following local variables are in scope:
#
#   ``PARENT_PACKAGE_NAME``
#
#     The name of the parent package.
#
#   ``SUBPACKAGE_NAME``
#
#     The local name of the subpackage (does not contain
#     the parent package name).
#
#   ``SUBPACKAGE_FULLNAME``
#
#     The full project-level name of the subpackage (which includes the parent
#     package name at the beginning,
#     ``${PARENT_PACKAGE_NAME}${SUBPACKAGE_NAME}``).
#
#   ``PACKAGE_NAME``
#
#     Inside the subpackage, the same as ``SUBPACKAGE_FULLNAME``.
#
macro(tribits_subpackage SUBPACKAGE_NAME_IN)

  if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("\nSUBPACKAGE: ${SUBPACKAGE_NAME_IN}")
  endif()

  tribits_subpackage_assert_call_context()

  # To provide context for various macros
  set(PACKAGE_NAME ${SUBPACKAGE_FULLNAME})

  set(PARENT_PACKAGE_SOURCE_DIR "${PACKAGE_SOURCE_DIR}")
  set(PARENT_PACKAGE_BINARY_DIR "${PACKAGE_BINARY_DIR}")

  # Now override the package-like variables
  tribits_set_common_vars(${SUBPACKAGE_FULLNAME})
  tribits_define_linkage_vars(${SUBPACKAGE_FULLNAME})
  tribits_pkg_init_exported_vars(${SUBPACKAGE_FULLNAME})

  tribits_append_package_specific_compiler_flags()
  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message("Final subpackage compiler flags:")
  endif()
  tribits_print_package_compiler_flags()

  # Set flags that are used  to check that macros are called in the correct order
  set(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED TRUE)
  set(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED FALSE)

endmacro()


function(tribits_subpackage_assert_call_context)

  # check that this is not being called from a package
  if (NOT CURRENTLY_PROCESSING_SUBPACKAGE)
    # we are in a package
    tribits_report_invalid_tribits_usage(
      "Cannot call tribits_subpackage() from a package."
      " Use tribits_package() instead"
      " ${CURRENT_PACKAGE_CMAKELIST_FILE}")
  else()
    # We are in a subpackage
    # check to see if postprocess is called before subpackage
    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "tribits_subpackage_postprocess() called before tribits_subpackage()")
    endif()

    # check to see if we have already called this macro
    if(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "Already called tribits_subpackge() for the"
        " ${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    endif()

    # make sure the name in the macro call matches the name in the packages cmake file
    if (NOT ${SUBPACKAGE_NAME_IN} STREQUAL ${SUBPACKAGE_NAME})
      tribits_report_invalid_tribits_usage(
        "Error, the package-defined subpackage name"
        " '${SUBPACKAGE_NAME_IN}' is not the same as the subpackage name"
        " '${SUBPACKAGE_NAME}' defined in the parent packages's"
        " Dependencies.cmake file")
    endif()
  endif()

endfunction()


# @MACRO: tribits_subpackage_postprocess()
#
# Macro that performs standard post-processing after defining a `TriBITS
# Subpackage`_ which is called at the bottom of a subpackage's
# `<packageDir>/<spkgDir>/CMakeLists.txt`_ file.
#
# Usage::
#
#   tribits_subpackage_postprocess()
#
# NOTE: This creates the aliased target ``${PACKAGE_NAME}::all_libs`` for all
# libraries in all subdirectories that don't have the TRIBITS_TESTONLY_LIB
# target property set on them.
#
# NOTE: It is unfortunate that a Subpackages's CMakeLists.txt file must call
# this macro but limitations of the CMake language make it necessary to do so.
#
macro(tribits_subpackage_postprocess)
  tribits_subpackage_postprocess_assert_call_context()
  tribits_package_postprocess_common()
endmacro()


macro(tribits_subpackage_postprocess_assert_call_context)

  # check that this is not being called from a package
  if (NOT CURRENTLY_PROCESSING_SUBPACKAGE)
    # This is being called from a package
    tribits_report_invalid_tribits_usage(
      "Cannot call tribits_subpackage_postprocess() from a package."
      " Use tribits_package_postprocess() instead"
      " ${CURRENT_PACKAGE_CMAKELIST_FILE}")
  else()
    # This is being called from a subpackage
    # check to make sure this has not already been called
    if (${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED)
      tribits_report_invalid_tribits_usage(
        "Already called tribits_subpackge_postprocess() for the"
        " ${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    endif()
  
    # make sure subpackage is called prior to subpackage postprocess
    if(NOT ${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_CALLED)
      tribits_report_invalid_tribits_usage(
        "tribits_subpackage() must be called before tribits_subpackage_postprocess()"
        " for the ${PARENT_PACKAGE_NAME} subpackage ${TRIBITS_SUBPACKAGE}")
    endif()
  endif()

  # Set flags that are used  to check that macros are called in the correct order
  dual_scope_set(${SUBPACKAGE_FULLNAME}_TRIBITS_SUBPACKAGE_POSTPROCESS_CALLED TRUE)

endmacro()
