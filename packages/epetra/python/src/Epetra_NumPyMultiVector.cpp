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

#include "Epetra_NumPyMultiVector.h"

#define DEBUG 0
#if DEBUG
#include <iostream>
#include <string>
using namespace std;
#endif

// Static variables
const Epetra_SerialComm   Epetra_NumPyMultiVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyMultiVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyMultiVector::tmp_map     = NULL               ;
PyArrayObject           * Epetra_NumPyMultiVector::tmp_range   = NULL               ;
char                    * Epetra_NumPyMultiVector::error_str   = NULL               ;

// Static helper functions
// =============================================================================
double * Epetra_NumPyMultiVector::getArray(PyObject * pyObject)
{
  // Try to build a contiguous PyArrayObject from the pyObject
  if (!tmp_array) tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject,'d',0,0);
  
  // If this fails, indicate an error by setting error_str, but build
  // a single vector with length zero to prevent a Bus Error
  if (!tmp_array) {
    error_str = "Error converting argument to an array";
    int dimensions[ ] = { 1, 0 };
    tmp_array = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);

  // If the contiguous PyArrayObject built successfully, make sure it has the correct
  // number of dimensions
  } else {
    if (tmp_array->nd < 2) {
      PyObject * tuple = Py_BuildValue("(ii)",1,-1);
      tmp_array = (PyArrayObject *) PyArray_Reshape(tmp_array,tuple);
      Py_DECREF(tuple);
    }
  }

  return (double *) (tmp_array->data);
}

// =============================================================================
double * Epetra_NumPyMultiVector::getArray(const Epetra_BlockMap & blockMap,
					   PyObject * pyObject)
{
  // PyObject argument is an int
  if PyInt_Check(pyObject) {
    int numVectors = (int) PyInt_AsLong(pyObject);
    int dimensions[ ] = { numVectors, blockMap.NumMyPoints() };
    tmp_array = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);
    if (tmp_array == NULL) {
      error_str  = "Error creating array";
      dimensions[0] = 1;
      dimensions[1] = 0;
      tmp_array  = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);
    }

  // PyObject argument is not an int ... try to build a contiguous PyArrayObject from it
  } else {
    if (!tmp_array) tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject,'d',0,0);

    // If this fails, indicate an error by setting error_str, but build
    // a single vector with length zero to prevent a Bus Error
    if (!tmp_array) {
      error_str = "Error converting argument to an array";
      int dimensions[ ] = { 1, 0 };
      tmp_array = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);

    // If the contiguous PyArrayObject built successfully, make sure it has the correct
    // number of dimensions
    } else {
      int   nd            = tmp_array->nd;
      int   dimensions[ ] = { 1, blockMap.NumMyPoints() };  // Default dimensions
      bool  reallocate    = false;
      int   arraySize     = _PyArray_multiply_list(tmp_array->dimensions,nd);

      if (nd < 2) {
	reallocate = true;
      } else {
	arraySize /= tmp_array->dimensions[0];
	if (arraySize != dimensions[1]) {
	  dimensions[0] = tmp_array->dimensions[0];
	  reallocate = true;
	}
      }

      // Reallocate the tmp_array if necessary
      if (reallocate) {
	PyArrayObject * myArray = (PyArrayObject *) PyArray_FromDims(2,dimensions,
								     PyArray_DOUBLE);
	double        * myData  = (double *) myArray->data;
	double        * tmpData = (double *) tmp_array->data;
	for (int i=0; i<dimensions[0]; i++) {
	  for (int j=0; j<arraySize && j<dimensions[1]; j++) {
	    myData[i*dimensions[1]+j] = tmpData[i*arraySize+j];
	  }
	}
	Py_XDECREF(tmp_array);
	tmp_array = myArray;
      }
    }
  }
  return (double *) (tmp_array->data);
}

// =============================================================================
Epetra_Map & Epetra_NumPyMultiVector::getMap(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  if (!tmp_map) {
    const int totalLength = PyArray_Size((PyObject *)tmp_array);
    tmp_map = new Epetra_Map(totalLength,0,defaultComm);
  }
  return *tmp_map;
}

// =============================================================================
int Epetra_NumPyMultiVector::getNumVectors(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return tmp_array->dimensions[0];
}

// =============================================================================
int Epetra_NumPyMultiVector::getNumVectors(const Epetra_BlockMap & blockMap,
					   PyObject * pyObject)
{
  if (!tmp_array) getArray(blockMap,pyObject);
  return tmp_array->dimensions[0];
}

// =============================================================================
int * Epetra_NumPyMultiVector::getRange(PyObject * range)
{
  // Try to create a contiguous array of integers from the PyObject
  if (!tmp_range) tmp_range = (PyArrayObject *) PyArray_ContiguousFromObject(range,'i',1,1);

  // If this fails, generate an array of integers of the same size as the PyObject
  if (tmp_range == NULL) {
    int dims[ ] = { PyObject_Length(range) };
    tmp_range = (PyArrayObject *) PyArray_FromDims(1,dims,PyArray_INT);
    int * data = (int *) tmp_range->data;
    for (int i=0; i<dims[0]; i++) data[i] = i;
  }

  // Obtain the length and return the array of integers
  return (int *) (tmp_range->data);
}

// =============================================================================
int Epetra_NumPyMultiVector::getRangeLen(PyObject * range)
{
  if (!tmp_range) getRange(range);
  return PyArray_Size((PyObject*)tmp_range);
}

// =============================================================================
int Epetra_NumPyMultiVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_map) getMap(pyObject);
  return tmp_map->NumMyPoints();
}

// =============================================================================

// Constructors
// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap,
						 int numVectors, bool zeroOut):
  Epetra_MultiVector(blockMap, numVectors, zeroOut)
{
  // Create the array object
  int dims[ ] = { numVectors, blockMap.NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(2,dims,PyArray_DOUBLE,
						    (char *)v[0]);

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);

  // Error message
  error_msg = error_str;
  error_str = NULL;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_MultiVector & source):
  Epetra_MultiVector(source)
{
  map = new Epetra_BlockMap(source.Map());
  int dims[ ] = { NumVectors(), map->NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(2,dims,PyArray_DOUBLE,
						    (char *)v[0]);

  // Error message
  error_msg = error_str;
  error_str = NULL;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap,
						 PyObject * pyObject):
  Epetra_MultiVector(View, blockMap, getArray(blockMap,pyObject),
		     blockMap.NumMyPoints(), getNumVectors(blockMap,pyObject))
{
  // Store the array pointer
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);

  // Error message
  error_msg = error_str;
  error_str = NULL;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(Epetra_DataAccess CV,
						 const Epetra_NumPyMultiVector & source,
						 PyObject * range):
  Epetra_MultiVector(CV, (const Epetra_MultiVector &) source, getRange(range),
		     getRangeLen(range))
{
  // Store the local map
  map = new Epetra_BlockMap(source.Map());

  // Inintialize the local Numeric array
  PyArrayObject * src_array = (PyArrayObject *) (source.ExtractView());
  int nd;
  // This shouldn't happen, but it does . . .
  if (NULL == src_array) nd = 2;
  else nd = src_array->nd;
  int dims[nd];
  dims[0] = NumVectors();
  if (NULL == src_array) dims[1] = source.MyLength();
  else for (int i=1; i<nd; i++) dims[i] = src_array->dimensions[i];

  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(nd,dims,PyArray_DOUBLE,
						    (char *)v[0]);

  // We're done with the tmp_range array
  Py_XDECREF(tmp_range);
  tmp_range = NULL;

  // Error message
  error_msg = error_str;
  error_str = NULL;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(PyObject * pyObject):
  Epetra_MultiVector(View, getMap(pyObject), getArray(pyObject),
		     getVectorSize(pyObject), getNumVectors(pyObject))
{
  // Store the pointer to the Epetra_Map
  map     = tmp_map;
  tmp_map = NULL;

  // Store the pointer to the PyArrayObject
  array     = tmp_array;
  tmp_array = NULL;

  // Error message
  error_msg = error_str;
  error_str = NULL;
}

// =============================================================================
// Destructor
Epetra_NumPyMultiVector::~Epetra_NumPyMultiVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::ErrorMsg() const
{
  if (error_msg == NULL) return Py_BuildValue("");
  return PyString_FromString(error_msg);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::ExtractCopy() const
{
  return PyArray_Copy(array);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Dot(const Epetra_MultiVector & a) const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::Dot(a, result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Dot returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm1() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;

  int status = Epetra_MultiVector::Norm1(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm1 returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*)(((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm2() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;

  int status = Epetra_MultiVector::Norm2(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm2 returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*)(((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::NormInf() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::NormInf(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "NormInf returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::NormWeighted(const Epetra_MultiVector & weights) const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::NormWeighted(weights,result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "NormWeighted returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MinValue() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::MinValue(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "MinValue returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MaxValue() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::MaxValue(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "MaxValue returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MeanValue() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};
  PyObject * po;
  double   * data;  

  int status = Epetra_MultiVector::MeanValue(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "MeanValue returned error code %d", status);
    goto fail;
  }
  po   = PyArray_FromDims(1, numVectors, PyArray_DOUBLE);
  data = (double*) (((PyArrayObject*)po)->data);
  for (int i=0; i<n; i++) data[i] = result[i];
  return po;
 fail:
  return NULL;
}

