# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

tribits_package_define_dependencies(
  LIB_OPTIONAL_PACKAGES Teuchos Epetra Zoltan Galeri Amesos Ifpack AztecOO EpetraExt Isorropia
  LIB_REQUIRED_TPLS BLAS LAPACK
  LIB_OPTIONAL_TPLS MPI METIS ParMETIS PETSC
  TEST_OPTIONAL_TPLS METIS ParMETIS
  )
