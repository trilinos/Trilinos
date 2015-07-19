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

#include "Epetra_NumPySerialSymDenseMatrix.hpp"

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
	PyArray_SimpleNew(2,dimensions,NPY_DOUBLE);
    }
    // If pyObject is not a bool nor an int, try to build a
    // contiguous 2D PyArrayObject from the pyObject
    else
    {
      // This function returns a borrowed ptr: no DECREF
      PyArray_Descr * dtype = 
	PyArray_DescrFromType(NPY_DOUBLE);
      tmp_array = (PyArrayObject *) PyArray_FromAny(pyObject, dtype, 2, 2,
						    NPY_ARRAY_FARRAY, NULL);
    }
  }
  // If no array has been correctly constructed, clean up and throw a
  // PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }

  return (double*)(PyArray_DATA(tmp_array));
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
    PyArray_Descr * dtype  = PyArray_DescrFromType(NPY_DOUBLE);
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
      double * oldData = Epetra_SerialSymDenseMatrix::A();
      double * newData = (double*) PyArray_DATA(array);
      int      size    = dimensions[0] * dimensions[1];
      for (int i=0; i<size; ++i) newData[i] = oldData[i];
    }
  }
}

// =============================================================================
int Epetra_NumPySerialSymDenseMatrix::getNumRows(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return PyArray_DIMS(tmp_array)[0];
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
