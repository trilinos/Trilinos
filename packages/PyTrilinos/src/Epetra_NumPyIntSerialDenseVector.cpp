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

#include "Epetra_NumPyIntSerialDenseVector.h"

namespace PyTrilinos
{

// Static variables
// =============================================================================
PyArrayObject * Epetra_NumPyIntSerialDenseVector::tmp_array = NULL;

// Static helper functions
// =============================================================================
int * Epetra_NumPyIntSerialDenseVector::getArray(PyObject * pyObject)
{
  // If tmp_array is NULL, build a PyArrayObject from the pyObject
  if (tmp_array == NULL)
  {
    // If pyObject is an int, build an array of that length
    if (PyInt_Check(pyObject))
    {
      npy_intp dimensions[ ] = {(npy_intp) PyInt_AsLong(pyObject)};
      tmp_array = (PyArrayObject*) PyArray_SimpleNew(1,dimensions,PyArray_INT);
    }
    // Else try to build a contiguous PyArrayObject from the pyObject
    else
    {
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,PyArray_INT,0,0);
    }
  }

  // If this fails, clean up and throw a PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }
  return (int*)(tmp_array->data);
}

// =============================================================================
void Epetra_NumPyIntSerialDenseVector::setArray()
{
  if (tmp_array)
  {
    array     = tmp_array;
    tmp_array = NULL;
  }
  else
  {
    npy_intp dimensions[ ] = {Length()};
    array = (PyArrayObject*)
      PyArray_SimpleNewFromData(1,dimensions,PyArray_INT,
				(void*)Epetra_IntSerialDenseVector::Values());
    if (!array)
    {
      cleanup();
      throw PythonException();
    }
  }
}

// =============================================================================
int Epetra_NumPyIntSerialDenseVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return (int) PyArray_MultiplyList(tmp_array->dimensions, tmp_array->nd);
}

void Epetra_NumPyIntSerialDenseVector::cleanup()
{
  if (tmp_array)
  {
    Py_DECREF(tmp_array);
    tmp_array = NULL;
  }
}

// Constructors
// =============================================================================
Epetra_NumPyIntSerialDenseVector::Epetra_NumPyIntSerialDenseVector():
  Epetra_IntSerialDenseVector()
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseVector::Epetra_NumPyIntSerialDenseVector(int length):
  Epetra_IntSerialDenseVector(length)
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseVector::Epetra_NumPyIntSerialDenseVector(PyObject * pyObject):
  Epetra_IntSerialDenseVector(View, getArray(pyObject), getVectorSize(pyObject))
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseVector
  setArray();
}

// =============================================================================
Epetra_NumPyIntSerialDenseVector::Epetra_NumPyIntSerialDenseVector(const Epetra_IntSerialDenseVector & src):
  Epetra_IntSerialDenseVector(src)
{
  // Synchronize the PyArrayObject with the Epetra_IntSerialDenseVector
  setArray();
}

// Destructor
// =============================================================================
Epetra_NumPyIntSerialDenseVector::~Epetra_NumPyIntSerialDenseVector()
{
  Py_XDECREF(array);
}

// Methods
// =============================================================================
int Epetra_NumPyIntSerialDenseVector::Size(int length)
{
  // Call the base-class method
  int result = Epetra_IntSerialDenseVector::Size(length);
  if (result)
  {
    PyErr_Format(PyExc_RuntimeError, "Size() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_IntSerialDenseVector data
  }
  return result;
}

// =============================================================================
int Epetra_NumPyIntSerialDenseVector::Resize(int length)
{
  // Call the base-class method
  int result = Epetra_IntSerialDenseVector::Resize(length);
  if (result)
  {
    PyErr_Format(PyExc_RuntimeError, "Resize() method failed with code %d", result);
  }
  else
  {
    Py_DECREF(array);   // Decrement the refcount to the current array
    setArray();         // Set the array from the Epetra_IntSerialDenseVector data
  }
  return result;
}

// =============================================================================
PyObject * Epetra_NumPyIntSerialDenseVector::Values() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

}  // Namespace PyTrilinos
