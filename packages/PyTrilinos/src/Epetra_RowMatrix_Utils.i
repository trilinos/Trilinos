// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%include "numpy.i"

%{
#include "Python.h"
#include "Epetra_RowMatrix.h"

// Helper functions for Epetra_RowMatrix
PyObject * Epetra_RowMatrix_GetEntries(const Epetra_RowMatrix& Matrix,
				       int GlobalRow) {
  int        myRow        = 0;
  int        numEntries   = 0;
  npy_intp   dims[1]      = { 0 };
  PyObject * indicesArray = NULL;
  PyObject * valuesArray  = NULL;
  int      * indices      = NULL;
  double   * values       = NULL;
  int        ierr         = 0;

  // Require Matrix be FillCompleted
  if (!Matrix.Filled()) {
    PyErr_SetString(PyExc_RuntimeError, "Matrix not FillCompleted");
    goto fail;
  }
  // Obtain the local row index from the global row index
  myRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  if (Matrix.NumMyRowEntries(myRow, numEntries)) {
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
  if (ierr != 0) {
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
				    int GlobalRow, int GlobalCol) {
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
  if (!Matrix.Filled()) {
    PyErr_SetString(PyExc_RuntimeError, "Matrix not FillCompleted");
    goto fail;
  }
  // Obtain the local row index from the global row index
  myRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  if (Matrix.NumMyRowEntries(myRow, numEntries)) {
    PyErr_Format(PyExc_ValueError, "Illegal global row index: %d", GlobalRow);
    goto fail;
  }
  // Obtain the local col index from the global col index
  myCol = Matrix.RowMatrixColMap().LID(GlobalCol);
  if (myCol < 0) {
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
  if (ierr != 0) {
    PyErr_Format(PyExc_RuntimeError, "ExtractMyRowCopy() returned %d error code", ierr);
    goto fail;
  }
  // Search for the requested Matrix entry
  for (int i = 0 ; i < numEntries ; ++i) {
    if (indices[i] == myCol) {
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
