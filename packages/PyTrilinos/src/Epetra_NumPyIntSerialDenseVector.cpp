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

#include "Epetra_NumPyIntSerialDenseVector.hpp"

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
      tmp_array = (PyArrayObject*) PyArray_SimpleNew(1,dimensions,NPY_INT);
    }
    // Else try to build a contiguous PyArrayObject from the pyObject
    else
    {
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,NPY_INT,0,0);
    }
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
      PyArray_SimpleNewFromData(1,dimensions,NPY_INT,
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
  return (int) PyArray_MultiplyList(PyArray_DIMS(tmp_array),
                                    PyArray_NDIM(tmp_array));
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
