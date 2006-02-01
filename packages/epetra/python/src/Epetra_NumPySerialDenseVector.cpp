// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Epetra_NumPySerialDenseVector.h"

// Static variable
// =============================================================================
PyArrayObject * Epetra_NumPySerialDenseVector::tmp_array = NULL;

// Static helper functions
// =============================================================================
double * Epetra_NumPySerialDenseVector::getArray(PyObject * pyObject)
{
  // If tmp_array is NULL, build a PyArrayObject from the pyObject
  if (!tmp_array) {
    // If pyObject is an int, build an array of that length
    if PyInt_Check(pyObject) {
      int dimensions[ ] = {(int) PyInt_AsLong(pyObject)};
      tmp_array = (PyArrayObject*) PyArray_FromDims(1,dimensions,'d');
    // Else try to build a contiguous PyArrayObject from the pyObject
    } else {
      tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject,'d',0,0);
    }
  }
  // If this fails, build a single vector with all zeros
  if (!tmp_array) {
    int dimensions[ ] = { PyObject_Length(pyObject) };
    tmp_array = (PyArrayObject *) PyArray_FromDims(1,dimensions,PyArray_DOUBLE);
  }

  return (double*)(tmp_array->data);
}

// =============================================================================
void Epetra_NumPySerialDenseVector::setArray()
{
  if (tmp_array) {
    array     = tmp_array;
    tmp_array = NULL;
  } else {
    int dimensions[ ] = {Length()};
    array = (PyArrayObject*) PyArray_FromDimsAndData(1,dimensions,'d',
						     (char*)Epetra_SerialDenseVector::Values());
  }
}

// =============================================================================
int Epetra_NumPySerialDenseVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return _PyArray_multiply_list(tmp_array->dimensions, tmp_array->nd);
}

// Constructors
// =============================================================================
Epetra_NumPySerialDenseVector::Epetra_NumPySerialDenseVector():
  Epetra_SerialDenseVector()
{
  // Synchronize the PyArrayObject with the Epetra_SerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPySerialDenseVector::Epetra_NumPySerialDenseVector(int length):
  Epetra_SerialDenseVector(length)
{
  // Synchronize the PyArrayObject with the Epetra_SerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPySerialDenseVector::Epetra_NumPySerialDenseVector(PyObject * pyObject):
  Epetra_SerialDenseVector(View, getArray(pyObject), getVectorSize(pyObject))
{
  // Synchronize the PyArrayObject with the Epetra_SerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPySerialDenseVector::Epetra_NumPySerialDenseVector(const Epetra_SerialDenseVector & src):
  Epetra_SerialDenseVector(src)
{
  // Synchronize the PyArrayObject with the Epetra_SerialDenseVector
  setArray();
}

// Destructor
// =============================================================================
Epetra_NumPySerialDenseVector::~Epetra_NumPySerialDenseVector() {
  Py_XDECREF(array);
}

// Methods
// =============================================================================
int Epetra_NumPySerialDenseVector::Size(int length) {
  // Call the base-class method
  int result = Epetra_SerialDenseVector::Size(length);
  if (result) {
    PyErr_Format(PyExc_RuntimeError, "Size() method failed with code %d", result);
  } else {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_SerialDenseVector data
  }
  return result;
}

// =============================================================================
int Epetra_NumPySerialDenseVector::Resize(int length) {
  // Call the base-class method
  int result = Epetra_SerialDenseVector::Resize(length);
  if (result) {
    PyErr_Format(PyExc_RuntimeError, "Resize() method failed with code %d", result);
  } else {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_SerialDenseVector data
  }
  return result;
}

// =============================================================================
PyObject * Epetra_NumPySerialDenseVector::Values() const {
  Py_INCREF(array);
  return PyArray_Return(array);
}
