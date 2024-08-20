# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

tribits_tpl_find_include_dirs_and_libraries( Pnetcdf
  REQUIRED_HEADERS pnetcdf.h
  REQUIRED_LIBS_NAMES pnetcdf
  )
