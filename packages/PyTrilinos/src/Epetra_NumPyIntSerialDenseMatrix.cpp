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

#include "Epetra_NumPyIntSerialDenseMatrix.hpp"

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
