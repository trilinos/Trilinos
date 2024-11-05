# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


# Only for unit testing
tribits_repository_define_tpls(
  TekoDepTPL    cmake/    PT
  )

# Above, the file cmake/FindTPLTekoDepTPL.cmake does not actually exit but it
# shows where a real extra repo might really put this file.
