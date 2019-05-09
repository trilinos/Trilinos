# -*- cmake -*-

# @HEADER
# ***********************************************************************
#
#          PyTrilinos: Python Interfaces to Trilinos Packages
#                 Copyright (2014) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact William F. Spotz (wfspotz@sandia.gov)
#
# ***********************************************************************
# @HEADER

ENABLE_TESTING()

INCLUDE(PyTrilinosMakeTest)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/testutil.py
  COMMAND ${PYTHON_EXECUTABLE}
          ${PyTrilinos_BINARY_DIR}/util/copyWithCMakeSubstitutions.py
          ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in
	  ${CMAKE_CURRENT_BINARY_DIR}/testutil.py
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in)

ADD_CUSTOM_TARGET(PyTrilinos_Example_testutil ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/testutil.py)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua
  COMMAND ${PYTHON_EXECUTABLE}
          ${PyTrilinos_BINARY_DIR}/util/copyWithCMakeSubstitutions.py
          ${CMAKE_CURRENT_SOURCE_DIR}/gre__115.mx
	  ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gre__115.mx)

ADD_CUSTOM_TARGET(PyTrilinos_gre__115 ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/gre__115.rua)

IF(PyTrilinos_ENABLE_Epetra)
  PyTrilinos_MAKE_TEST(exMPI_Init)
ENDIF(PyTrilinos_ENABLE_Epetra)

PyTrilinos_MAKE_TEST(exTeuchos)

IF(PyTrilinos_ENABLE_Epetra)
  PyTrilinos_MAKE_TEST(exEpetra                    )
  PyTrilinos_MAKE_TEST(exEpetra_Comm               )
  PyTrilinos_MAKE_TEST(exEpetra_CrsMatrix_Easy     )
  PyTrilinos_MAKE_TEST(exEpetra_CrsMatrix_Efficient)
  PyTrilinos_MAKE_TEST(exEpetra_FECrsMatrix_Easy   )
  PyTrilinos_MAKE_TEST(exEpetra_ImportExport       )
  IF(PyTrilinos_ENABLE_Tpetra)
    IF(PyTrilinos_ENABLE_Domi)
      PyTrilinos_MAKE_TEST(exEpetra_DAP)
    ENDIF(PyTrilinos_ENABLE_Domi)
  ENDIF(PyTrilinos_ENABLE_Tpetra)
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

IF(PyTrilinos_ENABLE_Domi)
  PyTrilinos_MAKE_TEST(exDomi_Simple)
ENDIF(PyTrilinos_ENABLE_Domi)

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
    IF(PyTrilinos_ENABLE_EpetraExt)
      IF(PyTrilinos_ENABLE_Triutils)
	IF(PyTrilinos_ENABLE_Galeri)
	  IF(PyTrilinos_ENABLE_AztecOO)
	    PyTrilinos_MAKE_TEST(exIFPACK)
	  ENDIF(PyTrilinos_ENABLE_AztecOO)
	ENDIF(PyTrilinos_ENABLE_Galeri)
      ENDIF(PyTrilinos_ENABLE_Triutils)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
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
    IF(PyTrilinos_ENABLE_EpetraExt)
      PyTrilinos_MAKE_TEST(exMLAPI_PyMatrix)
      PyTrilinos_MAKE_TEST(exMLAPI_Smoother)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
    IF(PyTrilinos_ENABLE_AztecOO)
      IF(NOT TPL_ENABLE_MPI)
	PyTrilinos_MAKE_TEST(exMLAPI_AztecOO)
      ENDIF(NOT TPL_ENABLE_MPI)
      IF(PyTrilinos_ENABLE_EpetraExt)
	IF(PyTrilinos_ENABLE_Galeri)
	  PyTrilinos_MAKE_TEST(exMLAPI_Simple)
	ENDIF(PyTrilinos_ENABLE_Galeri)
      ENDIF(PyTrilinos_ENABLE_EpetraExt)
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
      PyTrilinos_MAKE_TEST(exNOX_1Dfdm)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
  ENDIF(PyTrilinos_ENABLE_Epetra)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(NOX_ENABLE_LOCA)
  IF(PyTrilinos_ENABLE_NOX)
    IF(PyTrilinos_ENABLE_EpetraExt)
      IF(PyTrilinos_ENABLE_Epetra)
	PyTrilinos_MAKE_TEST(exLOCA_Chan     NUM_MPI_PROCS 1)
	#PyTrilinos_MAKE_TEST(exLOCA_Chan_Jac NUM_MPI_PROCS 1)
      ENDIF(PyTrilinos_ENABLE_Epetra)
    ENDIF(PyTrilinos_ENABLE_EpetraExt)
  ENDIF(PyTrilinos_ENABLE_NOX)
ENDIF(NOX_ENABLE_LOCA)
