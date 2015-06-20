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

#include "Epetra_NumPyIntVector.hpp"

#if NDARRAY_VERSION == 0x00090504
#define PyArray_ANYORDER -1
#endif

namespace PyTrilinos
{

// Static variables
const Epetra_SerialComm   Epetra_NumPyIntVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyIntVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyIntVector::tmp_map     = NULL               ;

// Static helper functions
// =============================================================================
int * Epetra_NumPyIntVector::getArray(PyObject * pyObject)
{
  // Try to build a contiguous PyArrayObject from the pyObject
  if (!tmp_array)
    tmp_array = (PyArrayObject *)
      PyArray_ContiguousFromObject(pyObject,NPY_INT,0,0);
  
  // If this fails, clean up and throw a PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }
  return (int*)(PyArray_DATA(tmp_array));
}

// =============================================================================
int * Epetra_NumPyIntVector::getArray(const Epetra_BlockMap & blockMap,
				      PyObject * pyObject)
{
  // Only build the tmp_array if it does not already exist
  if (!tmp_array)
  {
    // Default dimensions
    npy_intp defaultDims[ ] = { blockMap.NumMyPoints() };

    // PyObject argument is a bool
    if (PyBool_Check(pyObject))
    {
      tmp_array = (PyArrayObject *)
	PyArray_SimpleNew(1,defaultDims,NPY_INT);
    }
    // PyObject argument is not a bool ... try to build a contiguous
    // PyArrayObject from it
    else
    {
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,NPY_INT,0,0);
    }
    // If any PyArray factory functions fail, clean up and throw a
    // PythonException
    if (!tmp_array)
    {
      cleanup();
      throw PythonException();
    }
    int  nd = PyArray_NDIM(tmp_array);
    npy_intp arraySize = PyArray_MultiplyList(PyArray_DIMS(tmp_array),nd);
    if (arraySize != defaultDims[0])
    {
      PyArrayObject * myArray = (PyArrayObject *)
	PyArray_SimpleNew(1,defaultDims,NPY_INT);
      if (!myArray)
      {
	cleanup();
	throw PythonException();
      }
      int * myData  = (int *) PyArray_DATA(myArray);
      int * tmpData = (int *) PyArray_DATA(tmp_array);
      for (int i=0; i<defaultDims[0]; i++)
      {
	myData[i] = tmpData[i];
      }
      Py_XDECREF(tmp_array);
      tmp_array = myArray;
    }
  }
  return (int*)(PyArray_DATA(tmp_array));
}

// =============================================================================
Epetra_Map & Epetra_NumPyIntVector::getMap(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  if (!tmp_map)
  {
    const int totalLength = PyArray_Size((PyObject *)tmp_array);
    tmp_map = new Epetra_Map(totalLength,0,defaultComm);
  }
  return *tmp_map;
}

// =============================================================================
void Epetra_NumPyIntVector::cleanup()
{
  if (tmp_array)
  {
    Py_DECREF(tmp_array);
    tmp_array = NULL;
  }
  if (tmp_map)
  {
    delete tmp_map;
    tmp_map = NULL;
  }
}

// =============================================================================

// Constructors
// =============================================================================
Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_BlockMap & blockMap,
					     bool zeroOut):
  Epetra_IntVector(blockMap, zeroOut)
{
  // Create the array object
  npy_intp dims[ ] = { blockMap.NumMyPoints() };
  int *v = NULL;
  Epetra_IntVector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,NPY_INT,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_IntVector & source):
  Epetra_IntVector(source)
{
  map = new Epetra_BlockMap(source.Map());
  npy_intp dims[ ] = { map->NumMyPoints() };
  int *v = NULL;
  Epetra_IntVector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,NPY_INT,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
}

// =============================================================================
Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_BlockMap & blockMap,
					     PyObject * pyObject):
  Epetra_IntVector(View, blockMap, getArray(blockMap,pyObject))
{
  // Get the pointer to the array from static variable and clear
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyIntVector::Epetra_NumPyIntVector(PyObject * pyObject):
  Epetra_IntVector(View, getMap(pyObject), getArray(pyObject)) 
{
  // Store the pointer to the Epetra_Map
  assert(NULL != tmp_map);
  map     = tmp_map;
  tmp_map = NULL;

  // Store the pointer to the PyArrayObject
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;
}

// =============================================================================
// Destructor
Epetra_NumPyIntVector::~Epetra_NumPyIntVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyIntVector::ExtractCopy() const
{
  return PyArray_NewCopy(array,NPY_ANYORDER);
}

// =============================================================================
PyObject * Epetra_NumPyIntVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
PyObject * Epetra_NumPyIntVector::Values() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

}  // Namepace PyTrilinos
