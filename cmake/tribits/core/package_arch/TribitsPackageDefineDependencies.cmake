# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(TribitsGeneralMacros)


# @MACRO: tribits_package_define_dependencies()
#
# Define the dependencies for a given `TriBITS Package`_ (i.e. a top-level
# `TriBITS Package`_ or a `TriBITS Subpackage`_) in the package's
# `<packageDir>/cmake/Dependencies.cmake`_ file.
#
# Usage::
#
#   tribits_package_define_dependencies(
#      [LIB_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
#      [LIB_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
#      [TEST_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
#      [TEST_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
#      [LIB_REQUIRED_TPLS <tpl1> <tpl2> ...]
#      [LIB_OPTIONAL_TPLS <tpl1> <tpl2> ...]
#      [TEST_REQUIRED_TPLS <tpl1> <tpl2> ...]
#      [TEST_OPTIONAL_TPLS <tpl1> <tpl2> ...]
#      [SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
#        <spkg1_name>  <spkg1_dir>  <spkg1_classifications>  <spkg1_optreq>
#        <spkg2_name>  <spkg2_dir>  <spkg2_classifications>  <spkg2_optreq>
#        ...
#        ]
#      [REGRESSION_EMAIL_LIST  <regression-email-address>]
#      )
#
# Every argument in this macro is optional (that is, an package can have no
# upstream dependencies).  The arguments that apply to all packages are:
#
#   ``LIB_REQUIRED_PACKAGES``
#
#     List of required upstream packages that must be enabled in order to
#     build and use the libraries (or capabilities) in this package.
#
#   ``LIB_OPTIONAL_PACKAGES``
#
#     List of additional optional upstream packages that can be used in
#     this package if enabled.  These upstream packages need not be
#     enabled in order to enable this package but not enabling one or more
#     of these optional upstream packages will result in diminished
#     capabilities of this package.
#
#   ``TEST_REQUIRED_PACKAGES``
#
#     List of additional upstream packages that must be enabled in order to
#     build and/or run the tests and/or examples in this package.  If any of
#     these upstream packages are not enabled, then there will be no tests or
#     examples defined or run for this package.
#
#   ``TEST_OPTIONAL_PACKAGES``
#
#     List of additional optional upstream packages that can be used by the
#     tests in this package.  These upstream packages need not be enabled in
#     order to run some basic tests or examples for this package.  Typically,
#     extra tests that depend on optional test packages involve integration
#     testing of some type.  Not enabling these optional upstream packages
#     will result in diminished tests or examples.
#
#   ``LIB_REQUIRED_TPLS``
#
#     **DEPRECATED:** List of required upstream TPLs that must be enabled in
#     order to build and use the libraries (or capabilities) in this package.
#     (Add these to ``LIB_REQUIRED_PACKAGES`` instead.)
#
#   ``LIB_OPTIONAL_TPLS``
#
#     **DEPRECATED:** List of additional optional upstream TPLs that can be
#     used in this package if enabled.  These upstream TPLs need not be
#     enabled in order to use this package but not enabling one or more of
#     these optional upstream TPLs will result in diminished capabilities of
#     this package. (Add these to ``LIB_OPTIONAL_PACKAGES`` instead.)
#
#   ``TEST_REQUIRED_TPLS``
#
#     **DEPRECATED:** List of additional upstream TPLs that must be enabled in
#     order to build and/or run the tests and/or examples in this package.  If
#     any of these upstream TPLs are not enabled, then there will be no tests
#     or examples defined or run for this package.  (Add these to
#     ``TEST_REQUIRED_PACKAGES`` instead.)
#
#   ``TEST_OPTIONAL_TPLS``
#
#     **DEPRECATED:** List of additional optional upstream TPLs that can be
#     used by the tests in this package.  These upstream TPLs need not be
#     enabled in order to run basic tests for this package.  Typically, extra
#     tests that depend on optional TPLs involve integration testing or some
#     additional testing of some type.  (Add these to
#     ``TEST_OPTIONAL_PACKAGES`` instead.)
#
# NOTE: The above ``XXX_TPLS`` arguments/lists are **deprecated**.  At the
# package level, there is no distinction between upstream internal and
# external packages/TPLs, so all upstream package dependencies can (and should)
# be listed in the ``XXX_PACKAGES`` arguments/lists.  (There is no change in
# behavior listing upstream packages in ``XXX_PACKAGES`` or the ``XXX_TPLS``
# arguments/lists.)
#
# Only upstream packages can be listed (as defined by the order the packages
# are listed in `tribits_repository_define_packages()`_ in the
# `<repoDir>/PackagesList.cmake`_ or `<repoDir>/TPLsList.cmake`_ files).
# Otherwise an error will occur and processing will stop.  Misspelled package
# names are caught as well.
#
# Only direct package dependencies need to be listed.  Indirect package
# dependencies are automatically handled.  For example, if this package
# directly depends on package ``PKG2`` which depends on package ``PKG1``
# (but this package does not directly depend on anything in ``PKG1``) then
# this package only needs to list a dependency on ``PKG2``, not ``PKG1``.
# The dependency on ``PKG1`` will be taken care of automatically by the
# TriBITS dependency management system.
#
# The packages listed in ``LIB_REQUIRED_PACKAGES`` are implicitly also
# dependencies in ``TEST_REQUIRED_PACKAGES``.  Likewise
# ``LIB_OPTIONAL_PACKAGES`` are implicitly also dependencies in
# ``TEST_OPTIONAL_PACKAGES``.
#
# The upstream dependencies within a single list do not need to be listed in
# any particular order.  For example, if ``PKG2`` depends on ``PKG1``, and
# this given package depends on both, then one can list the dependencies as::
#
#   LIB_REQUIRED_PACKAGES PKG2 PKG1
#
# or::
#
#   LIB_REQUIRED_PACKAGES PKG1 PKG2
#
# If some upstream packages are allowed to be missing, this can be specified
# by calling the macro `tribits_allow_missing_external_packages()`_.
#
# A top-level `TriBITS Package`_ can also be broken down into `TriBITS
# Subpackages`_.  In this case, the following argument must be passed in:
#
#   .. _SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS:
#
#   ``SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS``
#
#     2D array with rows listing the subpackages where each row has the
#     columns:
#
#     * **SUBPACKAGE** (Column 0): The name of the subpackage ``<spkg_name>``.
#       The full package name is ``${PARENT_PACKAGE_NAME}<spkg_name>``.
#       The full package name is what is used in listing dependencies in
#       other packages.
#
#     * **DIRS** (Column 1): The subdirectory ``<spkg_dir>`` relative to the
#       parent package's base directory.  All of the contents of the
#       subpackage should be under this subdirectory.  This is assumed by the
#       TriBITS testing support software when mapping modified files to
#       packages that need to be tested (see `checkin-test.py`_).
#
#     * **CLASSIFICATIONS** (Column 2): The `Test Test Category`_ `PT`_,
#       `ST`_, `EX`_ and the maturity level ``EP``, ``RS``, ``PG``, ``PM``,
#       ``GRS``, ``GPG``, ``GPM``, and ``UM``, separated by a coma ',' with no
#       spaces in between (e.g. ``"PT,GPM"``).  These have exactly the same
#       meaning as for full packages (see
#       `tribits_repository_define_packages()`_).
#
#     * **OPTREQ** (Column 3): Determines if the outer parent package has an
#       ``OPTIONAL`` or ``REQUIRED`` dependence on this subpackage.
#
# Other variables that this macro handles:
#
#   ``REGRESSION_EMAIL_LIST``
#
#     The email list that is used to send CDash error messages.  If this
#     argument is missing, then the email list that CDash errors go to is
#     determined by other means (see `CDash regression email addresses`_).
#
macro(tribits_package_define_dependencies)

  cmake_parse_arguments(
     #prefix
     PARSE
     #options
     ""
     #one_value_keywords
     ""
     #multi_value_keywords
     "LIB_REQUIRED_PACKAGES;LIB_OPTIONAL_PACKAGES;TEST_REQUIRED_PACKAGES;TEST_OPTIONAL_PACKAGES;LIB_REQUIRED_TPLS;LIB_OPTIONAL_TPLS;TEST_REQUIRED_TPLS;TEST_OPTIONAL_TPLS;REGRESSION_EMAIL_LIST;SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS"
     ${ARGN}
     )

  tribits_check_for_unparsed_arguments()

  set(LIB_REQUIRED_DEP_PACKAGES ${PARSE_LIB_REQUIRED_PACKAGES})
  set(LIB_OPTIONAL_DEP_PACKAGES ${PARSE_LIB_OPTIONAL_PACKAGES})
  set(TEST_REQUIRED_DEP_PACKAGES ${PARSE_TEST_REQUIRED_PACKAGES})
  set(TEST_OPTIONAL_DEP_PACKAGES ${PARSE_TEST_OPTIONAL_PACKAGES})
  set(LIB_REQUIRED_DEP_TPLS ${PARSE_LIB_REQUIRED_TPLS})
  set(LIB_OPTIONAL_DEP_TPLS ${PARSE_LIB_OPTIONAL_TPLS})
  set(TEST_REQUIRED_DEP_TPLS ${PARSE_TEST_REQUIRED_TPLS})
  set(TEST_OPTIONAL_DEP_TPLS ${PARSE_TEST_OPTIONAL_TPLS})
  set(REGRESSION_EMAIL_LIST ${PARSE_REGRESSION_EMAIL_LIST})
  set(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    ${PARSE_SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS})

  # ToDo:
  # * Assert that REGRESSION_EMAIL_LIST has only one entry
  # * Assert that SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS is divisible
  #   by the number of columns!

endmacro()
