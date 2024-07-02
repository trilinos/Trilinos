# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Teuchos Teko
  LIB_OPTIONAL_PACKAGES Sacado
  LIB_REQUIRED_TPLS ExtraTPL2
  LIB_OPTIONAL_TPLS ExtraTPL1
  )
