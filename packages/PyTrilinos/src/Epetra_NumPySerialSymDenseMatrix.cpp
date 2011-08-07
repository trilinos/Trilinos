// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2010) Sandia Corporation
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

#include "Epetra_NumPySerialSymDenseMatrix.h"

namespace PyTrilinos
{

// Static variables
// =============================================================================
PyArrayObject * Epetra_NumPySerialSymDenseMatrix::tmp_array = NULL;

// Static helper functions
// =============================================================================
double * Epetra_NumPySerialSymDenseMatrix::getArray(PyObject * pyObject)
{
  // If tmp_array is NULL, build a PyArrayObject from the pyObject
  if (!tmp_array)
  {
    // If pyObject is an int, then emulate an int-int constructor
    if (PyInt_Check(pyObject))
    {
      int numRows = (int) PyInt_AsLong(pyObject);
      npy_intp dimensions[ ] = {numRows, numRows};
      tmp_array = (PyArrayObject *)
	PyArray_SimpleNew(2,dimensions,PyArray_DOUBLE);
    }
    // If pyObject is not a bool nor an int, try to build a
    // contiguous 2D PyArrayObject from the pyObject
    else
    {
      // This function returns a borrowed ptr: no DECREF
      PyArray_Descr * dtype = 
	PyArray_DescrFromType(PyArray_DOUBLE);
      tmp_array = (PyArrayObject *) PyArray_FromAny(pyObject, dtype, 2, 2,
						    NPY_FARRAY, NULL);
    }
  }
  // If no array has been correctly constructed, clean up and throw a
  // PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }

  return (double*)(tmp_array->data);
}

// =============================================================================
void Epetra_NumPySerialSymDenseMatrix::setArray(bool copy)
{
  if (tmp_array)
  {
    array     = tmp_array;
    tmp_array = NULL;
  }
  else
  {
    npy_intp dimensions[ ] = { M(), N() };
    double * data = NULL;
    if (!copy) data = Epetra_SerialSymDenseMatrix::A();
    // This NumPy function returns a borrowed pointer: no DECREF
    PyArray_Descr * dtype  = PyArray_DescrFromType(PyArray_DOUBLE);
    array = (PyArrayObject*)
      PyArray_NewFromDescr(&PyArray_Type,dtype,2,dimensions,NULL,(void*)data,
			   NPY_FARRAY,NULL);
    if (!array)
    {
      cleanup();
      throw PythonException();
    }
    if (copy)
    {
      double * oldData = Epetra_SerialSymDenseMatrix::A();
      double * newData = (double*) array->data;
      int      size    = dimensions[0] * dimensions[1];
      for (int i=0; i<size; ++i) newData[i] = oldData[i];
    }
  }
}

// =============================================================================
int Epetra_NumPySerialSymDenseMatrix::getNumRows(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return tmp_array->dimensions[0];
}

// =============================================================================
void Epetra_NumPySerialSymDenseMatrix::cleanup()
{
  if (tmp_array)
  {
    Py_DECREF(tmp_array);
    tmp_array = NULL;
  }
}

// Constructors
// =============================================================================
Epetra_NumPySerialSymDenseMatrix::Epetra_NumPySerialSymDenseMatrix() :
  Epetra_SerialSymDenseMatrix()
{
  // Synchronize the PyArrayObject with the Epetra_SerialSymDenseMatrix
  setArray();
}

// =============================================================================
Epetra_NumPySerialSymDenseMatrix::Epetra_NumPySerialSymDenseMatrix(PyObject * pyObject) :
  Epetra_SerialSymDenseMatrix(View,
			      getArray(  pyObject),
			      getNumRows(pyObject),
			      getNumRows(pyObject) )
{
  // Synchronize the PyArrayObject with the Epetra_SerialSymDenseMatrix
  setArray();
}

// =============================================================================
Epetra_NumPySerialSymDenseMatrix::Epetra_NumPySerialSymDenseMatrix(
  const Epetra_SerialSymDenseMatrix & src) :
  Epetra_SerialSymDenseMatrix(src)
{
  // Synchronize the PyArrayObject with the Epetra_SerialSymDenseMatrix
  setArray(true);
}

// Destructor
// =============================================================================
Epetra_NumPySerialSymDenseMatrix::~Epetra_NumPySerialSymDenseMatrix()
{
  Py_XDECREF(array);
}

// Methods
// =============================================================================
double Epetra_NumPySerialSymDenseMatrix::operator() (int rowIndex, int colIndex)
{
  return double(Epetra_SerialSymDenseMatrix::operator()(rowIndex, colIndex));
}

// =============================================================================
int Epetra_NumPySerialSymDenseMatrix::Shape(int numRowsCols)
{
  // Call the base-class method
  int result = Epetra_SerialSymDenseMatrix::Shape(numRowsCols);
  if (result)
  {
    PyErr_Format(PyExc_ValueError, "Shape() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_SerialSymDenseMatrix data
  }
  return result;
}

// =============================================================================
int Epetra_NumPySerialSymDenseMatrix::Reshape(int numRowsCols)
{
  // Call the base-class method
  int result = Epetra_SerialSymDenseMatrix::Reshape(numRowsCols);
  if (result)
  {
    PyErr_Format(PyExc_ValueError, "Reshape() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_SerialSymDenseMatrix data
  }
  return result;
}

// =============================================================================
PyObject * Epetra_NumPySerialSymDenseMatrix::A()
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

}  // Namespace PyTrilinos
