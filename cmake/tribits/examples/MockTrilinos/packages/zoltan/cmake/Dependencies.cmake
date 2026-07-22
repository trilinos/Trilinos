# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

if (SHOW_INVALID_TPL_NAME_ERROR)
  set(Zlib_str Zlib)
else()
  set(Zlib_str "")
endif()

tribits_package_define_dependencies(
  LIB_OPTIONAL_TPLS MPI ParMETIS Scotch ${Zlib_str}
  )
