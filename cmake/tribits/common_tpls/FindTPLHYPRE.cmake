# @HEADER
# *****************************************************************************
#           Trilinos: An Object-Oriented Solver Framework
#
# Copyright 2001-2024 NTESS and the Trilinos contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER


TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HYPRE
  REQUIRED_HEADERS HYPRE.h HYPRE_config.h
  REQUIRED_LIBS_NAMES HYPRE
  )
