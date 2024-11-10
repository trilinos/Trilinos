# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

if (EXTRA_REPO_INCLUDE_MISSING_OPTIONAL_DEP_PACKAGE)
   set(MISSING_OPTIONAL_DEP_PACKAGE MissingUpstreamPackage)
  tribits_allow_missing_external_packages(MissingUpstreamPackage)
else()
   set(MISSING_OPTIONAL_DEP_PACKAGE)
endif()

if (EXTRA_REPO_INCLUDE_MISSING_REQUIRED_DEP_PACKAGE)
   set(MISSING_REQUIRED_DEP_PACKAGE MissingUpstreamPackage)
  tribits_allow_missing_external_packages(MissingUpstreamPackage)
else()
   set(MISSING_REQUIRED_DEP_PACKAGE)
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Teuchos ${MISSING_REQUIRED_DEP_PACKAGE}
  LIB_OPTIONAL_PACKAGES ${MISSING_OPTIONAL_DEP_PACKAGE}
  LIB_REQUIRED_TPLS Boost
  REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov
  )


