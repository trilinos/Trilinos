# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

if (EXTRA_REPO_PACKAGE2_ADD_TEST_DEPS)
  set(ExRepo2Package2_TEST_DEPS
    TEST_REQUIRED_PACKAGES Teuchos
    TEST_OPTIONAL_PACKAGES RTOp Ex2Package1
    TEST_REQUIRED_TPLS Boost
    TEST_OPTIONAL_TPLS MPI Boost
    )
  # NOTE: Include duplicates above so that we can test removing them as well!
else()
  set(ExRepo2Package2_TEST_DEPS "")
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Teuchos
  LIB_OPTIONAL_PACKAGES Ex2Package1
  ${ExRepo2Package2_TEST_DEPS}
  )
