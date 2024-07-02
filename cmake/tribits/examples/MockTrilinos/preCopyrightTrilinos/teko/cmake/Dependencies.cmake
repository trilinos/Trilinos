# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

set(INVALID_TESTING_PACKAGE_NAME)
if (SHOW_INVALID_PACKAGE_NAME_ERROR)
  set(INVALID_TESTING_PACKAGE_NAME Anasazi)
endif()

set(SELF_REFERENCE_TESTING_PACKAGE_NAME)
if (SHOW_SELF_PACKAGE_DEPENDENCY_ERROR)
  set(SELF_REFERENCE_TESTING_PACKAGE_NAME Teko)
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Teuchos Epetra Thyra EpetraExt Stratimikos
    AztecOO ${INVALID_TESTING_PACKAGE_NAME} ML
    ${SELF_REFERENCE_TESTING_PACKAGE_NAME} Ifpack Amesos
  LIB_OPTIONAL_PACKAGES Isorropia
  LIB_REQUIRED_TPLS TekoDepTPL
  )
