# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


macro(tribits_repository_setup_extra_options)

  # Needed for testing
  set(TPL_ENABLE_MPI OFF CACHE BOOL "Enable MPI support.")

  # Just an example of a repo-wide option
  advanced_set(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )

endmacro()
