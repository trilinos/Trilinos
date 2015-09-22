// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%include "numpy.i"

// Helper functions here use some code in numpy.i fragments, so we
// nedd to activate those fragments.
%fragment("NumPy_Fragments");

%{
#include "Python.h"
#include "Epetra_RowMatrix.h"

// Helper functions for Epetra_RowMatrix
PyObject * Epetra_RowMatrix_GetEntries(const Epetra_RowMatrix& Matrix,
				       int GlobalRow)
{
  int        myRow        = 0;
  int        numEntries   = 0;
  npy_intp   dims[1]      = { 0 };
  PyObject * indicesArray = NULL;
  PyObject * valuesArray  = NULL;
  int      * indices      = NULL;
  double   * values       = NULL;
  int        ierr         = 0;

  // Require Matrix be FillCompleted
  if (!Matrix.Filled())
  {
    PyErr_SetString(PyExc_RuntimeError, "Matrix not FillCompleted");
    goto fail;
  }
  // Obtain the local row index from the global row index
  myRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  if (Matrix.NumMyRowEntries(myRow, numEntries))
  {
    PyErr_Format(PyExc_ValueError, "Illegal global row index: %d", GlobalRow);
    goto fail;
  }
  // Build the NumPy arrays for the indices and values
  dims[0] = numEntries;
  indicesArray = PyArray_SimpleNew(1, dims, NPY_INT   );
  if (indicesArray == NULL) goto fail;
  valuesArray  = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  if (valuesArray == NULL) goto fail;
  indices = (int   *) array_data(indicesArray);
  values  = (double*) array_data(valuesArray );
  //Extract the row data
  ierr = Matrix.ExtractMyRowCopy(myRow, (int)dims[0], numEntries, values, indices);
  if (ierr != 0)
  {
    PyErr_Format(PyExc_RuntimeError, "ExtractMyRowCopy() returned %d error code", ierr);
    goto fail;
  }
  // Convert the local row indexes to global row indexes
  for (int i = 0 ; i < numEntries ; ++i)
    indices[i] = Matrix.RowMatrixColMap().GID(indices[i]);
  // Return the arrays
  return Py_BuildValue("(OO)", indicesArray, valuesArray);
 fail:
  Py_XDECREF(indicesArray);
  Py_XDECREF(valuesArray);
  return NULL;
}

PyObject* Epetra_RowMatrix_GetEntry(const Epetra_RowMatrix& Matrix,
				    int GlobalRow, int GlobalCol)
{
  int        myRow        = 0;
  int        myCol        = 0;
  int        numEntries   = 0;
  npy_intp   dims[1]      = { 0 };
  PyObject * indicesArray = NULL;
  PyObject * valuesArray  = NULL;
  int      * indices      = NULL;
  double   * values       = NULL;
  int        ierr         = 0;
  double     value        = 0.0;

  // Require Matrix be FillCompleted
  if (!Matrix.Filled())
  {
    PyErr_SetString(PyExc_RuntimeError, "Matrix not FillCompleted");
    goto fail;
  }
  // Obtain the local row index from the global row index
  myRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  if (Matrix.NumMyRowEntries(myRow, numEntries))
  {
    PyErr_Format(PyExc_ValueError, "Illegal global row index: %d", GlobalRow);
    goto fail;
  }
  // Obtain the local col index from the global col index
  myCol = Matrix.RowMatrixColMap().LID(GlobalCol);
  if (myCol < 0)
  {
    PyErr_Format(PyExc_ValueError, "Illegal global col index: %d", GlobalCol);
    goto fail;
  }
  // Build the NumPy arrays for the indices and values
  dims[0] = numEntries;
  indicesArray = PyArray_SimpleNew(1, dims, NPY_INT   );
  if (indicesArray == NULL) goto fail;
  valuesArray  = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  if (valuesArray == NULL) goto fail;
  indices = (int   *) array_data(indicesArray);
  values  = (double*) array_data(valuesArray );
  //Extract the row data
  ierr = Matrix.ExtractMyRowCopy(myRow, (int)dims[0], numEntries, values, indices);
  if (ierr != 0)
  {
    PyErr_Format(PyExc_RuntimeError, "ExtractMyRowCopy() returned %d error code", ierr);
    goto fail;
  }
  // Search for the requested Matrix entry
  for (int i = 0 ; i < numEntries ; ++i)
  {
    if (indices[i] == myCol)
    {
      value = values[i];
      break;
    }
  }
  // Clean up the arrays and return the result
  Py_DECREF(indicesArray);
  Py_DECREF(valuesArray );
  return PyFloat_FromDouble(value);
 fail:
  Py_XDECREF(indicesArray);
  Py_XDECREF(valuesArray);
  return NULL;
}

%}
