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
  COMMAND ${PyTrilinos_BINARY_DIR}/util/copyWithCMakeSubstitutions.py
          ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in
	  ${CMAKE_CURRENT_BINARY_DIR}/testutil.py
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/testutil.py.in)

ADD_CUSTOM_TARGET(PyTrilinos_Test_testutil ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/testutil.py)

PyTrilinos_MAKE_TEST(testTeuchos              )
PyTrilinos_MAKE_TEST(testTeuchos_Comm         )
PyTrilinos_MAKE_TEST(testTeuchos_ParameterList)
PyTrilinos_MAKE_TEST(testTeuchos_RCP          )
PyTrilinos_MAKE_TEST(testTeuchos_ScalarTraits )

IF(PyTrilinos_ENABLE_Epetra)
  PyTrilinos_MAKE_TEST(testEpetra_BlockMap      )
  PyTrilinos_MAKE_TEST(testEpetra_Comm          )
  PyTrilinos_MAKE_TEST(testEpetra_CompObject    )
  PyTrilinos_MAKE_TEST(testEpetra_CrsGraph      )
  PyTrilinos_MAKE_TEST(testEpetra_CrsMatrix     )
  PyTrilinos_MAKE_TEST(testEpetra_FEVector      )
  PyTrilinos_MAKE_TEST(testEpetra_ImportExport  )
  PyTrilinos_MAKE_TEST(testEpetra_IntSerialDense)
  PyTrilinos_MAKE_TEST(testEpetra_IntVector     )
  PyTrilinos_MAKE_TEST(testEpetra_LinearProblem )
  PyTrilinos_MAKE_TEST(testEpetra_LocalMap      )
  PyTrilinos_MAKE_TEST(testEpetra_Map           )
  PyTrilinos_MAKE_TEST(testEpetra_MapColoring   )
  PyTrilinos_MAKE_TEST(testEpetra_MultiVector   )
  PyTrilinos_MAKE_TEST(testEpetra_Object        )
  PyTrilinos_MAKE_TEST(testEpetra_SerialDense   )
  PyTrilinos_MAKE_TEST(testEpetra_Time          )
  PyTrilinos_MAKE_TEST(testEpetra_Util          )
  PyTrilinos_MAKE_TEST(testEpetra_Vector        )
ENDIF(PyTrilinos_ENABLE_Epetra)

IF(PyTrilinos_ENABLE_Triutils)
  PyTrilinos_MAKE_TEST(testTriUtils)
ENDIF(PyTrilinos_ENABLE_Triutils)

IF(PyTrilinos_ENABLE_Tpetra)
  PyTrilinos_MAKE_TEST(testTpetra_Map        )
  PyTrilinos_MAKE_TEST(testTpetra_Vector     )
  PyTrilinos_MAKE_TEST(testTpetra_MultiVector)
ENDIF(PyTrilinos_ENABLE_Tpetra)

IF(PyTrilinos_ENABLE_EpetraExt)
  PyTrilinos_MAKE_TEST(testEpetraExt_Coloring)
  PyTrilinos_MAKE_TEST(testEpetraExt_Property)
  IF(EpetraExt_USING_HDF5)
    PyTrilinos_MAKE_TEST(testEpetraExt_HDF5)
  ENDIF(EpetraExt_USING_HDF5)
ENDIF(PyTrilinos_ENABLE_EpetraExt)

IF(PyTrilinos_ENABLE_Isorropia)
  PyTrilinos_MAKE_TEST(testIsorropia)
ENDIF(PyTrilinos_ENABLE_Isorropia)

IF(PyTrilinos_ENABLE_AztecOO)
  PyTrilinos_MAKE_TEST(testAztecOO)
ENDIF(PyTrilinos_ENABLE_AztecOO)

IF(PyTrilinos_ENABLE_Galeri)
  PyTrilinos_MAKE_TEST(testGaleri_CrsMatrices)
  PyTrilinos_MAKE_TEST(testGaleri_Maps       )
ENDIF(PyTrilinos_ENABLE_Galeri)

IF(PyTrilinos_ENABLE_Amesos)
  PyTrilinos_MAKE_TEST(testAmesos)
ENDIF(PyTrilinos_ENABLE_Amesos)

IF(PyTrilinos_ENABLE_Ifpack)
  IF(PyTrilinos_ENABLE_EpetraExt)
    PyTrilinos_MAKE_TEST(testIFPACK)
  ENDIF(PyTrilinos_ENABLE_EpetraExt)
ENDIF(PyTrilinos_ENABLE_Ifpack)

IF(PyTrilinos_ENABLE_Komplex)
  PyTrilinos_MAKE_TEST(testKomplex_LinearProblem)
ENDIF(PyTrilinos_ENABLE_Komplex)

IF(PyTrilinos_ENABLE_Anasazi)
  PyTrilinos_MAKE_TEST(testAnasazi)
ENDIF(PyTrilinos_ENABLE_Anasazi)

IF(PyTrilinos_ENABLE_Pliris)
  PyTrilinos_MAKE_TEST(testPliris)
ENDIF(PyTrilinos_ENABLE_Pliris)

IF(PyTrilinos_ENABLE_ML)
  IF(PyTrilinos_ENABLE_EpetraExt)
    PyTrilinos_MAKE_TEST(testML_Space)
    PyTrilinos_MAKE_TEST(testML_CompObject)
    PyTrilinos_MAKE_TEST(testML_TimeObject)
    PyTrilinos_MAKE_TEST(testML_MultiVector)
    PyTrilinos_MAKE_TEST(testML_PyMatrix)
    PyTrilinos_MAKE_TEST(testML_Preconditioner)
  ENDIF(PyTrilinos_ENABLE_EpetraExt)
ENDIF(PyTrilinos_ENABLE_ML)

IF(PyTrilinos_ENABLE_NOX)
  PyTrilinos_MAKE_TEST(testNOX_StatusTest)
ENDIF(PyTrilinos_ENABLE_NOX)

IF(PyTrilinos_ENABLE_LOCA)
ENDIF(PyTrilinos_ENABLE_LOCA)
