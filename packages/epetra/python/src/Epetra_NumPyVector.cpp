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

#include "float.h"

#include "Epetra_NumPyVector.h"

#define DEBUG 0
#if DEBUG
#include <iostream>
#include <string>
using namespace std;
#endif

// Static variables
const Epetra_SerialComm   Epetra_NumPyVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyVector::tmp_map     = NULL               ;

// Static helper functions
// =============================================================================
Epetra_Map & Epetra_NumPyVector::getEpetraMap(PyObject * pyObject)
{
  getSourceData(pyObject);   // This creates the tmp_array
  const int totalLength = PyArray_Size((PyObject *)tmp_array);
  assert(NULL == tmp_map);
  tmp_map = new Epetra_Map(totalLength,0,defaultComm);
  return *tmp_map;
}

// =============================================================================
double * Epetra_NumPyVector::getSourceData(PyObject * pyObject)
{
  assert(NULL == tmp_array);
  tmp_array = (PyArrayObject *) PyArray_ContiguousFromObject(pyObject, 'd', 0, 0);
  return (double *)tmp_array->data;
}
// =============================================================================

// Constructors
// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap & blockMap, bool zeroOut):
  Epetra_Vector(blockMap, zeroOut)
{
  // Create the array object
  int dims[ ] = { blockMap.NumMyElements() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,
						    (char *)v);

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_NumPyVector & source):
  Epetra_Vector(source)
{
  map = new Epetra_BlockMap(source.Map());
  int dims[ ] = { map->NumMyElements() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,
						    (char *)v);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap & blockMap,
                                       PyObject * pyObject):
  Epetra_Vector(View, blockMap, getSourceData(pyObject))
{
  // Get the pointer to the array from static variable and clear
  assert(NULL != tmp_array);
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
// Epetra_NumPyVector::Epetra_NumPyVector(Epetra_DataAccess CV, const Epetra_MultiVector & source,
// 				       int index):
//   Epetra_Vector(CV,source,index)
// {
//   map = new Epetra_BlockMap(source.Map());
//   int dims[ ] = { map->NumMyElements() };
//   double *v = NULL;
//   ExtractView(&v);
//   array = (PyArrayObject *) PyArray_FromDimsAndData(1,dims,PyArray_DOUBLE,
// 						    (char *)v);
// }

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(PyObject * pyObject):
  //Epetra_Vector(View, getEpetraMap(pyObject), (double*)tmp_array->data) 
  Epetra_Vector(View, getEpetraMap(pyObject), getSourceData(pyObject)) 
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
Epetra_NumPyVector::~Epetra_NumPyVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyVector::ExtractCopy() const
{
  return PyArray_Copy(array);
}

// =============================================================================
PyObject * Epetra_NumPyVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
PyObject * Epetra_NumPyVector::getArray()
{
  return PyArray_Return(array);
}

// =============================================================================
double Epetra_NumPyVector::Norm1() const {
  double result[1];
  int status = Epetra_MultiVector::Norm1(result);
  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm1 returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyVector::Norm2() const {
  double result[1];
  int status = Epetra_MultiVector::Norm2(result);
  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Norm2 returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyVector::NormInf() const {
  double result[1];
  int status = Epetra_MultiVector::NormInf(result);
  if (status) {
    PyErr_Format(PyExc_RuntimeError, "NormInf returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyVector::Dot(const Epetra_Vector & A) const {
  double result[1];
  int status = Epetra_MultiVector::Dot(A, result);
  if (status) {
    PyErr_Format(PyExc_RuntimeError, "Dot returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}
