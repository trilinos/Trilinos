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

#include "Epetra_NumPyIntSerialDenseMatrix.h"

namespace PyTrilinos
{

// Static variables
// =============================================================================
PyArrayObject * Epetra_NumPyIntSerialDenseMatrix::tmp_array = NULL;

// Static helper functions
// =============================================================================
int * Epetra_NumPyIntSerialDenseMatrix::getArray(PyObject * pyObject)
{
  // If tmp_array is NULL, build a two-dimensional PyArrayObject from the pyObject
  if (tmp_array == NULL)
  {
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_INT);
    tmp_array = (PyArrayObject *)
      PyArray_FromAny(pyObject,dtype,2,2,NPY_ARRAY_FARRAY,NULL);
  }

  // If this fails, clean up and throw a PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }

  return (int*)(PyArray_DATA(tmp_array));
}

// =============================================================================
void Epetra_NumPyIntSerialDenseMatrix::setArray(bool copy)
{
  if (tmp_array)
  {
    array     = tmp_array;
    tmp_array = NULL;
  }
  else
  {
    npy_intp dimensions[ ]  = { M(), N() };
    int    * data           = NULL;
    if (!copy) data = Epetra_IntSerialDenseMatrix::A();
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_INT);
    array = (PyArrayObject*)
      PyArray_NewFromDescr(&PyArray_Type,dtype,2,dimensions,NULL,(void*)data,
			   NPY_ARRAY_FARRAY,NULL);
    if (!array)
    {
      cleanup();
      throw PythonException();
    }
    if (copy)
    {
      int * oldData = Epetra_IntSerialDenseMatrix::A();
      int * newData = (int*) PyArray_DATA(array);
      int   size    = dimensions[0] * dimensions[1];
      for (int i=0; i<size; ++i) newData[i] = oldData[i];
    }
  }
}

// =============================================================================
int Epetra_NumPyIntSerialDenseMatrix::getNumCols(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return (int) PyArray_DIMS(tmp_array)[1];
}

// =============================================================================
int Epetra_NumPyIntSerialDenseMatrix::getNumRows(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return (int) PyArray_DIMS(tmp_array)[0];
}

// =============================================================================
void Epetra_NumPyIntSerialDenseMatrix::cleanup()
{
  if (tmp_array)
  {
    Py_DECREF(tmp_array);
    tmp_array = NULL;
  }
}

// Constructors
// =============================================================================
Epetra_NumPyIntSerialDenseMatrix::Epetra_NumPyIntSerialDenseMatrix():
  Epetra_IntSerialDenseMatrix()
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseMatrix
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseMatrix::Epetra_NumPyIntSerialDenseMatrix(int numRows, int numCols):
  Epetra_IntSerialDenseMatrix(numRows, numCols)
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseMatrix
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseMatrix::Epetra_NumPyIntSerialDenseMatrix(PyObject * pyObject):
  Epetra_IntSerialDenseMatrix(View, getArray(pyObject), getNumRows(pyObject),
			      getNumRows(pyObject), getNumCols(pyObject))
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseMatrix
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseMatrix::Epetra_NumPyIntSerialDenseMatrix(
  const Epetra_IntSerialDenseMatrix & src):
  Epetra_IntSerialDenseMatrix(src)
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseMatrix
  setArray(true);
}

// Destructor
// =============================================================================
Epetra_NumPyIntSerialDenseMatrix::~Epetra_NumPyIntSerialDenseMatrix()
{
  Py_XDECREF(array);
}

// Methods
// =============================================================================
int Epetra_NumPyIntSerialDenseMatrix::operator() (int rowIndex, int colIndex)
{
  return int(Epetra_IntSerialDenseMatrix::operator()(rowIndex,colIndex));
}

// =============================================================================
int Epetra_NumPyIntSerialDenseMatrix::Shape(int numRows, int numCols)
{
  // Call the base-class method
  int result = Epetra_IntSerialDenseMatrix::Shape(numRows, numCols);
  if (result)
  {
    PyErr_Format(PyExc_ValueError, "Shape() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_IntSerialDenseMatrix data
  }
  return result;
}

// =============================================================================
int Epetra_NumPyIntSerialDenseMatrix::Reshape(int numRows, int numCols)
{
  // Call the base-class method
  int result = Epetra_IntSerialDenseMatrix::Reshape(numRows, numCols);
  if (result)
  {
    PyErr_Format(PyExc_ValueError, "Reshape() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_IntSerialDenseMatrix data
  }
  return result;
}

// =============================================================================
PyObject * Epetra_NumPyIntSerialDenseMatrix::A()
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

}  // Namespace PyTrilinos
