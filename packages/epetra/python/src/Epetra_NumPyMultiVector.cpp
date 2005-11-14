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
const Epetra_SerialComm   Epetra_NumPyMultiVector::defaultComm 	 = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyMultiVector::tmp_array   	 = NULL               ;
Epetra_Map              * Epetra_NumPyMultiVector::tmp_map     	 = NULL               ;
PyArrayObject           * Epetra_NumPyMultiVector::tmp_range   	 = NULL               ;
int                       Epetra_NumPyMultiVector::tmp_range_len = 0                  ;

// Static helper functions
// =============================================================================
Epetra_Map & Epetra_NumPyMultiVector::getEpetraMap(PyObject * pyObject)
{
  getArrayFromObject(pyObject);   // This creates the tmp_array
  const int totalLength = PyArray_Size((PyObject *)tmp_array);
  assert(NULL == tmp_map);
  tmp_map = new Epetra_Map(totalLength,0,defaultComm);
  return *tmp_map;
}

// =============================================================================
int * Epetra_NumPyMultiVector::getRange(PyObject * range)
{
  assert(NULL == tmp_range    );
  assert(0    == tmp_range_len);

  // Try to create a contiguous array of integers from the PyObject
  tmp_range = (PyArrayObject *) PyArray_ContiguousFromObject(range,'i',1,1);

  // If this fails, generate an array of integers of the same size as the PyObject
  if (tmp_range == NULL) {
    int dims[ ] = { PyObject_Length(range) };
    tmp_range = (PyArrayObject *) PyArray_FromDims(1,dims,PyArray_INT);
    int * data = (int *) tmp_range->data;
    for (int i=0; i<dims[0]; i++) data[i] = i;
  }

  // Obtain the length and return the array of integers
  tmp_range_len = PyArray_Size((PyObject *)tmp_range);
  return (int *) (tmp_range->data);
}

// =============================================================================
double * Epetra_NumPyMultiVector::getArrayFromObject(PyObject * pyObject)
{
  assert(NULL == tmp_array);

  // Try to build a contiguous PyArrayObject from the pyObject
  tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject,'d',0,0);
  
  // If this fails, build a single vector with all zeros
  if (tmp_array == NULL) {
    int dimensions[ ] = { 1, PyObject_Length(pyObject) };
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
double * Epetra_NumPyMultiVector::getArrayFromMapAndObject(const Epetra_BlockMap & blockMap,
							   PyObject * pyObject)
{
  assert(NULL == tmp_array);

  // PyObject argument is an int
  if PyInt_Check(pyObject) {
    int numVectors = (int) PyInt_AsLong(pyObject);
    int dimensions[ ] = { numVectors, blockMap.NumMyElements() };
    tmp_array = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);

  // PyObject argument is not an int ... try to build a contiguous PyArrayObject from it
  } else {
    tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject,'d',0,0);

    // If this fails, build a single vector with all zeros
    if (tmp_array == NULL) {
      int dimensions[ ] = { 1, blockMap.NumMyElements() };
      tmp_array = (PyArrayObject *) PyArray_FromDims(2,dimensions,PyArray_DOUBLE);

    // If the contiguous PyArrayObject built successfully, make sure it has the correct
    // number of dimensions
    } else {
      int   nd            = tmp_array->nd;
      int   dimensions[ ] = { 1, blockMap.NumMyElements() };  // Default dimensions
      bool  reallocate    = false;
      int   arraySize     = 1;
      for (int i=0; i<tmp_array->nd; i++) arraySize *= tmp_array->dimensions[i];

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

// Constructors
// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap,
						 int numVectors, bool zeroOut):
  Epetra_MultiVector(blockMap, numVectors, zeroOut)
{
  // Create the array object
  int dims[ ] = { numVectors, blockMap.NumMyElements() };
  double **v = NULL;
  ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(2,dims,PyArray_DOUBLE,
						    (char *)v[0]);

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_NumPyMultiVector & source):
  Epetra_MultiVector(source)
{
  map = new Epetra_BlockMap(source.Map());
  int dims[ ] = { NumVectors(), map->NumMyElements() };
  double **v = NULL;
  ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(2,dims,PyArray_DOUBLE,
						    (char *)v[0]);
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap,
						 PyObject * pyObject):
  Epetra_MultiVector(View, blockMap, getArrayFromMapAndObject(blockMap,pyObject),
		     blockMap.NumMyElements(), tmp_array->dimensions[0])
{
  // Get the pointer to the array from static variable and clear
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);

  // Perform some checks
  double ** eptr;
  double * nptr;
  ExtractView(&eptr);
  nptr = (double *) (array->data);
  assert(eptr[0] == nptr);

  printf("Epetra  pointer is %d\n",int(eptr[0]));
  printf("Numeric pointer is %d\n",int(nptr   ));

}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(Epetra_DataAccess CV,
						 const Epetra_NumPyMultiVector & source,
						 PyObject * range):
  Epetra_MultiVector(CV, (const Epetra_MultiVector &) source, getRange(range),
		     tmp_range_len)
{
  // Store the local map
  map = new Epetra_BlockMap(source.Map());

  // Inintialize the local Numeric array
  PyArrayObject * src_array = (PyArrayObject *) (source.getArray());
  int nd;
  // This shouldn't happen, but it does . . .
  if (NULL == src_array) nd = 2;
  else nd = src_array->nd;
  int dims[nd];
  dims[0] = NumVectors();
  if (NULL == src_array) dims[1] = source.MyLength();
  else for (int i=1; i<nd; i++) dims[i] = src_array->dimensions[i];

  double **v = NULL;
  ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(nd,dims,PyArray_DOUBLE,
						    (char *)v[0]);

  // We're done with the tmp_range array
  assert(NULL != tmp_range);
  Py_XDECREF(tmp_range);
  tmp_range = NULL;
  assert(0 != tmp_range_len);
  tmp_range_len = 0;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(PyObject * pyObject):
  Epetra_MultiVector(View, getEpetraMap(pyObject), getArrayFromObject(pyObject),
		     tmp_map->NumMyElements(), tmp_array->dimensions[0]) 
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
Epetra_NumPyMultiVector::~Epetra_NumPyMultiVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::getArray() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm1() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};

  int status = Epetra_MultiVector::Norm1(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm1 returned error code %d", status);
    goto fail;
  }
  return PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE, (char *)result);
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm2() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};

  int status = Epetra_MultiVector::Norm2(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm2 returned error code %d", status);
    goto fail;
  }
  return PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE, (char *)result);
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::NormInf() const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};

  int status = Epetra_MultiVector::NormInf(result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "NormInf returned error code %d", status);
    goto fail;
  }
  return PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE, (char *)result);
 fail:
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Dot(const Epetra_MultiVector & A) const {
  int    n = NumVectors();
  double result[n];
  int    numVectors[ ] = {n};

  int status = Epetra_MultiVector::Dot(A, result);

  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Dot returned error code %d", status);
    goto fail;
  }
  return PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE, (char *)result);
 fail:
  return NULL;
}
