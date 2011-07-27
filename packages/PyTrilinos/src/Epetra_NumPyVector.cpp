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

#include "Epetra_NumPyVector.h"

#if NDARRAY_VERSION == 0x00090504
#define PyArray_ANYORDER -1
#endif

namespace PyTrilinos
{

// Static variables
const Epetra_SerialComm   Epetra_NumPyVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyVector::tmp_map     = NULL               ;

// Static helper functions
// =============================================================================
double * Epetra_NumPyVector::getArray(PyObject * pyObject)
{
  // Try to build a contiguous PyArrayObject from the pyObject
  if (!tmp_array)
    tmp_array = (PyArrayObject *)
      PyArray_ContiguousFromObject(pyObject,PyArray_DOUBLE,0,0);
  
  // If this fails, clean up and throw a PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }

  return (double*)(tmp_array->data);
}

// =============================================================================
double * Epetra_NumPyVector::getArray(const Epetra_BlockMap & blockMap,
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
	PyArray_SimpleNew(1,defaultDims,PyArray_DOUBLE);
      if (tmp_array == NULL)
      {
	cleanup();
	throw PythonException();
      }
      else
      {
	if (pyObject == Py_True)     // bool zeroOut is True
	{
	  double * data = (double*) tmp_array->data;
	  for (int i=0; i<defaultDims[0]; ++i) data[i] = 0.0;
	}
      }

    // PyObject argument is not an bool ... try to build a contiguous
    // PyArrayObject from it
    }
    else
    {
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,PyArray_DOUBLE,0,0);

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
	int nd = tmp_array->nd;
	npy_intp arraySize = PyArray_MultiplyList(tmp_array->dimensions,nd);
	if (arraySize != defaultDims[0])
	{
	  PyArrayObject * myArray = (PyArrayObject *)
	    PyArray_SimpleNew(1,defaultDims,PyArray_DOUBLE);
	  double * myData  = (double *) myArray->data;
	  double * tmpData = (double *) tmp_array->data;
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
  return (double*)(tmp_array->data);
}

// =============================================================================
Epetra_Map & Epetra_NumPyVector::getMap(PyObject * pyObject)
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
int Epetra_NumPyVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_map) getMap(pyObject);
  return tmp_map->NumMyPoints();
}

// =============================================================================
void Epetra_NumPyVector::cleanup()
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
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap & blockMap, bool zeroOut):
  Epetra_Vector(blockMap, zeroOut)
{
  // Create the array object
  npy_intp dims[ ] = { blockMap.NumMyPoints() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_Vector & source):
  Epetra_Vector(source)
{
  map = new Epetra_BlockMap(source.Map());
  npy_intp dims[ ] = { map->NumMyPoints() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap & blockMap,
                                       PyObject * pyObject):
  Epetra_Vector(View, blockMap, getArray(blockMap,pyObject))
{
  // Get the pointer to the array from static variable and clear
  array     = tmp_array;
  tmp_array = NULL;

  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(Epetra_DataAccess CV,
				       const Epetra_Vector & source):
  Epetra_Vector(CV,source,0)
{
  // Store the Epetra_Vector's map
  map = new Epetra_BlockMap(source.Map());

  // Wrap the Epetra_Vector
  npy_intp dims[ ] = { (npy_intp) map->NumMyElements() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(Epetra_DataAccess CV,
				       const Epetra_MultiVector & source,
				       int index):
  Epetra_Vector(CV,source,index)
{
  // Store the Epetra_MultiVector's map
  map = new Epetra_BlockMap(source.Map());

  // Wrap the Epetra_MultiVector
  npy_intp dims[ ] = { map->NumMyElements() };
  double *v = NULL;
  Epetra_Vector::ExtractView(&v);
  array = (PyArrayObject *)
    PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(void *)v);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
}

// =============================================================================
Epetra_NumPyVector::Epetra_NumPyVector(PyObject * pyObject):
  Epetra_Vector(View, getMap(pyObject), getArray(pyObject)) 
{
  // Store the pointer to the Epetra_Map
  map     = tmp_map;
  tmp_map = NULL;

  // Store the pointer to the PyArrayObject
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
  return PyArray_NewCopy(array,PyArray_ANYORDER);
}

// =============================================================================
PyObject * Epetra_NumPyVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
double Epetra_NumPyVector::Dot(const Epetra_Vector & A) const
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
double Epetra_NumPyVector::Norm1() const
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
double Epetra_NumPyVector::Norm2() const
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
double Epetra_NumPyVector::NormInf() const
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
double Epetra_NumPyVector::NormWeighted(const Epetra_Vector & weights) const
{
  double result[1];
  int status = Epetra_MultiVector::NormWeighted(weights,result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "NormWeighted returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
double Epetra_NumPyVector::MinValue() const
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
double Epetra_NumPyVector::MaxValue() const
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
double Epetra_NumPyVector::MeanValue() const
{
  double result[1];
  int status = Epetra_MultiVector::MeanValue(result);
  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MeanValue returned error code %d", status);
    goto fail;
  }
  return result[0];
 fail:
  // Return type is double, so we cannot return NULL on failure
  return DBL_MAX;
}

// =============================================================================
int Epetra_NumPyVector::ReplaceGlobalValues(PyObject * values, PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::ReplaceGlobalValues(lenValues,
					      (double *) myValues->data,
					      (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::ReplaceGlobalValues(int blockOffset, PyObject * values,
					    PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::ReplaceGlobalValues(lenValues, blockOffset,
					      (double *) myValues->data,
					      (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::ReplaceMyValues(PyObject * values, PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::ReplaceMyValues(lenValues,
					  (double *) myValues->data,
					  (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::ReplaceMyValues(int blockOffset, PyObject * values,
					PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::ReplaceMyValues(lenValues, blockOffset,
					  (double *) myValues->data,
					  (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::SumIntoGlobalValues(PyObject * values, PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::SumIntoGlobalValues(lenValues,
					      (double *) myValues->data,
					      (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::SumIntoGlobalValues(int blockOffset, PyObject * values,
					    PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::SumIntoGlobalValues(lenValues, blockOffset,
					      (double *) myValues->data,
					      (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::SumIntoMyValues(PyObject * values, PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::SumIntoMyValues(lenValues,
					  (double *) myValues->data,
					  (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

// =============================================================================
int Epetra_NumPyVector::SumIntoMyValues(int blockOffset, PyObject * values,
					PyObject * indices)
{
  PyArrayObject * myValues  = NULL;
  PyArrayObject * myIndices = NULL;
  int lenValues;
  int lenIndices;
  int result;
  myValues  = (PyArrayObject *)
    PyArray_ContiguousFromObject(values,PyArray_DOUBLE,0,0);
  if (!myValues) goto fail;
  myIndices = (PyArrayObject *)
    PyArray_ContiguousFromObject(indices,PyArray_INT,0,0);
  if (!myIndices) goto fail;
  lenValues  = (int) PyArray_MultiplyList(myValues->dimensions,  myValues->nd );
  lenIndices = (int) PyArray_MultiplyList(myIndices->dimensions, myIndices->nd);
  if (lenValues != lenIndices)
  {
    PyErr_Format(PyExc_ValueError,
		 "Sequence lengths are %d and %d, but must be of same length",
		 lenValues, lenIndices);
    goto fail;
  }
  result = Epetra_Vector::SumIntoMyValues(lenValues, blockOffset,
					  (double *) myValues->data,
					  (int    *) myIndices->data);
  Py_DECREF(myValues );
  Py_DECREF(myIndices);
  return result;
 fail:
  Py_XDECREF(myValues );
  Py_XDECREF(myIndices);
  return -1;
}

}  // Namespace PyTrilinos
