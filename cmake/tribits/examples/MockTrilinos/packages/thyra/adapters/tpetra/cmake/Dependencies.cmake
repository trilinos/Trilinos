# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

set(SELF_SUBPACKAGE_DEP_REFERENCE)
if (SHOW_SELF_SUBPACKAGE_DEPENDENCY_ERROR)
  set(SELF_SUBPACKAGE_DEP_REFERENCE ThyraTpetra)
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Tpetra ${SELF_SUBPACKAGE_DEP_REFERENCE} ThyraCoreLibs
  )
