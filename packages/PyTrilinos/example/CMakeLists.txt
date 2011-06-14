# -*- cmake -*-
# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2010) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

ENABLE_TESTING()

INCLUDE(PyTrilinosMakeTest)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/testutil.py
  COMMAND ${CMAKE_SOURCE_DIR}/packages/PyTrilinos/util/copyWithCMakeSubstitutions.py
          ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in
	  ${CMAKE_CURRENT_BINARY_DIR}/testutil.py
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in)

ADD_CUSTOM_TARGET(ex_testutil ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/testutil.py)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua
  COMMAND ${CMAKE_SOURCE_DIR}/packages/PyTrilinos/util/copyWithCMakeSubstitutions.py
          ${CMAKE_CURRENT_SOURCE_DIR}/gre__115.mx
	  ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gre__115.mx)

ADD_CUSTOM_TARGET(gre__115 ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua)

IF(PyTrilinos_ENABLE_Teuchos)
  IF(PyTrilinos_ENABLE_Epetra)
    PyTrilinos_MAKE_TEST(exMPI_Init)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Teuchos)
  PyTrilinos_MAKE_TEST(exTeuchos)
ENDIF(PyTrilinos_ENABLE_Teuchos)

IF(PyTrilinos_ENABLE_Epetra)
  PyTrilinos_MAKE_TEST(exEpetra                    )
  PyTrilinos_MAKE_TEST(exEpetra_Comm               )
  PyTrilinos_MAKE_TEST(exEpetra_CrsMatrix_Easy     )
  PyTrilinos_MAKE_TEST(exEpetra_CrsMatrix_Efficient)
  PyTrilinos_MAKE_TEST(exEpetra_FECrsMatrix_Easy   )
  PyTrilinos_MAKE_TEST(exEpetra_ImportExport       )
ENDIF(PyTrilinos_ENABLE_Epetra)

IF(PyTrilinos_ENABLE_TriUtils)
  PyTrilinos_MAKE_TEST(exTriUtils)
ENDIF(PyTrilinos_ENABLE_TriUtils)

IF(PyTrilinos_ENABLE_EpetraExt)
  IF(PyTrilinos_ENABLE_Epetra)
    PyTrilinos_MAKE_TEST(exEpetraExt_IO_MatrixMarket)
    PyTrilinos_MAKE_TEST(exEpetraExt_MatrixMatrix   )
    IF(EpetraExt_USING_HDF5)
      PyTrilinos_MAKE_TEST(exEpetraExt_IO_HDF5)
    ENDIF(EpetraExt_USING_HDF5)
    # Until expat is integrated as an official TPL, the following
    # example script is doomed to failure
    #PyTrilinos_MAKE_TEST(exEpetraExt_IO_XML)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_Isorropia)
  IF(PyTrilinos_ENABLE_Epetra)
    PyTrilinos_MAKE_TEST(exIsorropia)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Isorropia)

IF(PyTrilinos_ENABLE_AztecOO)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(PyTrilinos_ENABLE_Galeri)
      PyTrilinos_MAKE_TEST(exAztecOO)
    ENDIF(PyTrilinos_ENABLE_Galeri)
    PyTrilinos_MAKE_TEST(exAztecOO_Operator
      XHOST godel.sandia.gov)
    PyTrilinos_MAKE_TEST(exAztecOO_RowMatrix
      XHOST godel.sandia.gov)
    PyTrilinos_MAKE_TEST(exAztecOO_BasicRowMatrix)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  IF(PyTrilinos_ENABLE_Epetra)
    PyTrilinos_MAKE_TEST(exGaleri       )
    PyTrilinos_MAKE_TEST(exGaleri_ReadHB)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  IF(PyTrilinos_ENABLE_Epetra)
    PyTrilinos_MAKE_TEST(exAmesos_Factory)
    PyTrilinos_MAKE_TEST(exAmesos_Solvers)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(PyTrilinos_ENABLE_Triutils)
      IF(PyTrilinos_ENABLE_Galeri)
	IF(PyTrilinos_ENABLE_AztecOO)
	  PyTrilinos_MAKE_TEST(exIFPACK)
	ENDIF(PyTrilinos_ENABLE_AztecOO)
      ENDIF(PyTrilinos_ENABLE_Galeri)
    ENDIF(PyTrilinos_ENABLE_Triutils)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(PyTrilinos_ENABLE_AztecOO)
      PyTrilinos_MAKE_TEST(exKomplex_Simple)
    ENDIF(PyTrilinos_ENABLE_AztecOO)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(PyTrilinos_ENABLE_Galeri)
      PyTrilinos_MAKE_TEST(exAnasazi_BlockDavidson)
    ENDIF(PyTrilinos_ENABLE_Galeri)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(NOT TPL_ENABLE_MPI)
      PyTrilinos_MAKE_TEST(exMLAPI         )
      PyTrilinos_MAKE_TEST(exMLAPI_Iterate )
    ENDIF(NOT TPL_ENABLE_MPI)
    PyTrilinos_MAKE_TEST(exMLAPI_PyMatrix)
    PyTrilinos_MAKE_TEST(exMLAPI_Smoother)
    IF(PyTrilinos_ENABLE_AztecOO)
      IF(NOT TPL_ENABLE_MPI)
	PyTrilinos_MAKE_TEST(exMLAPI_AztecOO)
      ENDIF(NOT TPL_ENABLE_MPI)
      IF(PyTrilinos_ENABLE_Galeri)
	PyTrilinos_MAKE_TEST(exMLAPI_Simple)
      ENDIF(PyTrilinos_ENABLE_Galeri)
    ENDIF(PyTrilinos_ENABLE_AztecOO)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  IF(PyTrilinos_ENABLE_Epetra)
    IF(NOT TPL_ENABLE_MPI)
      PyTrilinos_MAKE_TEST(exNOX_2DSim)
      IF(PyTrilinos_ENABLE_EpetraExt)
	PyTrilinos_MAKE_TEST(exNOXEpetraExt_2DSim)
      ENDIF(PyTrilinos_ENABLE_EpetraExt)
    ENDIF(NOT TPL_ENABLE_MPI)
    IF(PyTrilinos_ENABLE_EpetraExt)
      IF(PyTrilinos_ENABLE_Teuchos)
	PyTrilinos_MAKE_TEST(exNOX_1Dfdm)
      ENDIF(PyTrilinos_ENABLE_Teuchos)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(PyTrilinos_ENABLE_LOCA)
  IF(PyTrilinos_ENABLE_NOX)
    IF(PyTrilinos_ENABLE_EpetraExt)
      IF(PyTrilinos_ENABLE_Epetra)
	IF(PyTrilinos_ENABLE_Teuchos)
	  #PyTrilinos_MAKE_TEST(exLOCA_Chan)
	ENDIF(PyTrilinos_ENABLE_Teuchos)
      ENDIF(PyTrilinos_ENABLE_Epetra)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
  ENDIF(PyTrilinos_ENABLE_NOX)
ENDIF(PyTrilinos_ENABLE_LOCA)
