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

#include "float.h"

#include "Epetra_NumPyFEVector.h"

#if NDARRAY_VERSION == 0x00090504
#define PyArray_ANYORDER -1
#endif

namespace PyTrilinos
{

// Static variables
const Epetra_SerialComm   Epetra_NumPyFEVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyFEVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyFEVector::tmp_map     = NULL               ;

// Static helper functions
// =============================================================================
double * Epetra_NumPyFEVector::getArray(PyObject * pyObject)
{
  // Try to build a contiguous PyArrayObject from the pyObject
  if (!tmp_array)
    tmp_array = (PyArrayObject *)
      PyArray_ContiguousFromObject(pyObject,NPY_DOUBLE,0,0);
  
  // If this fails, build a vector with length zero to prevent a Bus
  // Error
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }

  return (double*)(PyArray_DATA(tmp_array));
}

// =============================================================================
double * Epetra_NumPyFEVector::getArray(const Epetra_BlockMap & blockMap,
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
	PyArray_SimpleNew(1,defaultDims,NPY_DOUBLE);
      if (tmp_array == NULL)
      {
	cleanup();
	throw PythonException();
      }
      else
      {
	if (pyObject == Py_True)     // bool zeroOut is True
	{
	  double * data = (double*) PyArray_DATA(tmp_array);
	  for (int i=0; i<defaultDims[0]; ++i) data[i] = 0.0;
	}
      }
    }
    // PyObject argument is not an bool ... try to build a contiguous
    // PyArrayObject from it
    else
    {
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,NPY_DOUBLE,0,0);

      // If this fails, clean up and throw a PythonException
      if (!tmp_array)
      {
	cleanup();
	throw PythonException();
      }
      // If the contiguous PyArrayObject built successfully, make sure
      // it has the correct number of dimensions
      else
      {
	int  nd = PyArray_NDIM(tmp_array);
	npy_intp arraySize = PyArray_MultiplyList(PyArray_DIMS(tmp_array),nd);
	if (arraySize != defaultDims[0])
	{
	  PyArrayObject * myArray = 
	    (PyArrayObject *) PyArray_SimpleNew(1,defaultDims,NPY_DOUBLE);
	  if (!myArray)
	  {
	    cleanup();
	    throw PythonException();
	  }
	  double        * myData  = (double *) PyArray_DATA(myArray);
	  double        * tmpData = (double *) PyArray_DATA(tmp_array);
	  for (int i=0; i<defaultDims[0]; i++)
	  {
	    myData[i] = tmpData[i];
	  }
	  Py_XDECREF(tmp_array);
	  tmp_array = myArray;
	}
      }
    }
  }
  return (double*)(PyArray_DATA(tmp_array));
}

// =============================================================================
Epetra_Map & Epetra_NumPyFEVector::getMap(PyObject * pyObject)
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
int Epetra_NumPyFEVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_map) getMap(pyObject);
  return tmp_map->NumMyPoints();
}

// =============================================================================
void Epetra_NumPyFEVector::cleanup()
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
Epetra_NumPyFEVector::Epetra_NumPyFEVector(const Epetra_BlockMap & blockMap,
					   int numVectors,
					   bool ignoreNonLocalEntries):
  Epetra_FEVector(blockMap, numVectors, ignoreNonLocalEntries)
{
  // Create the array object
  npy_intp dims[ ] = { numVectors, blockMap.NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,(void *)v[0]);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyFEVector::Epetra_NumPyFEVector(const Epetra_FEVector & source):
  Epetra_FEVector(source)
{
  map = new Epetra_BlockMap(source.Map());
  npy_intp dims[ ] = { source.NumVectors(), map->NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,(void *)v[0]);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
}

// =============================================================================
// Destructor
Epetra_NumPyFEVector::~Epetra_NumPyFEVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyFEVector::ExtractCopy() const
{
  return PyArray_NewCopy(array,NPY_ANYORDER);
}

// =============================================================================
PyObject * Epetra_NumPyFEVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
double Epetra_NumPyFEVector::Dot(const Epetra_FEVector & A) const
{
  double result[1];
  int status = Epetra_MultiVector::Dot(A, result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Dot returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::Norm1() const
{
  double result[1];
  int status = Epetra_MultiVector::Norm1(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Norm1 returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::Norm2() const
{
  double result[1];
  int status = Epetra_MultiVector::Norm2(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Norm2 returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::NormInf() const
{
  double result[1];
  int status = Epetra_MultiVector::NormInf(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "NormInf returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::NormWeighted(const Epetra_FEVector & weights) const
{
  double result[1];
  int status = Epetra_MultiVector::NormWeighted(weights,result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "NormWeighted returned error code %d",
		 status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::MinValue() const
{
  double result[1];
  int status = Epetra_MultiVector::MinValue(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MinValue returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::MaxValue() const
{
  double result[1];
  int status = Epetra_MultiVector::MaxValue(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MaxValue returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyFEVector::MeanValue() const
{
  double result[1];
  int status = Epetra_MultiVector::MeanValue(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MeanValue returned error code %d",
		 status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
int Epetra_NumPyFEVector::ReplaceGlobalValues(PyObject * indices,
					      PyObject * values)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values, NPY_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,NPY_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(PyArray_DIMS(myValues),
                                          PyArray_NDIM(myValues) );
  lenIndices = (int) PyArray_MultiplyList(PyArray_DIMS(myIndices),
                                          PyArray_NDIM(myIndices));
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_FEVector::ReplaceGlobalValues(lenValues,
						(int    *) PyArray_DATA(myIndices),
						(double *) PyArray_DATA(myValues));
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyFEVector::SumIntoGlobalValues(PyObject * indices,
					      PyObject * values)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,NPY_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,NPY_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(PyArray_DIMS(myValues),
                                          PyArray_NDIM(myValues) );
  lenIndices = (int) PyArray_MultiplyList(PyArray_DIMS(myIndices),
                                          PyArray_NDIM(myIndices));
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_FEVector::SumIntoGlobalValues(lenValues,
						(int    *) PyArray_DATA(myIndices),
						(double *) PyArray_DATA(myValues));
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

}  // Namespace PyTrilinos
