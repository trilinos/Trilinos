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

#include "Epetra_NumPyMultiVector.h"

#if NDARRAY_VERSION == 0x00090504
#define PyArray_ANYORDER -1
#endif

namespace PyTrilinos
{

// Static variables
const Epetra_SerialComm   Epetra_NumPyMultiVector::defaultComm = Epetra_SerialComm();
PyArrayObject           * Epetra_NumPyMultiVector::tmp_array   = NULL               ;
Epetra_Map              * Epetra_NumPyMultiVector::tmp_map     = NULL               ;
PyArrayObject           * Epetra_NumPyMultiVector::tmp_range   = NULL               ;

// Static helper functions
// =============================================================================
double * Epetra_NumPyMultiVector::getArray(PyObject * pyObject)
{
  // Try to build a contiguous PyArrayObject from the pyObject
  if (!tmp_array)
    tmp_array = (PyArrayObject *)
      PyArray_ContiguousFromObject(pyObject,NPY_DOUBLE,0,0);
  
  // If this fails, clean up and throw a PythonException
  if (!tmp_array)
  {
    cleanup();
    throw PythonException();
  }
  // If the contiguous PyArrayObject built successfully, make sure it has the correct
  // number of dimensions
  else
  {
    if (PyArray_NDIM(tmp_array) < 2)
    {
      PyObject * tuple = Py_BuildValue("(ii)",1,-1);
      tmp_array = (PyArrayObject *) PyArray_Reshape(tmp_array,tuple);
      Py_DECREF(tuple);
    }
  }

  return (double *) (PyArray_DATA(tmp_array));
}

// =============================================================================
double * Epetra_NumPyMultiVector::getArray(const Epetra_BlockMap & blockMap,
					   PyObject * pyObject)
{
  // PyObject argument is an int
  if (PyInt_Check(pyObject))
  {
    int numVectors = (int) PyInt_AsLong(pyObject);
    npy_intp dimensions[ ] = { numVectors, blockMap.NumMyPoints() };
    tmp_array = (PyArrayObject *)
      PyArray_SimpleNew(2,dimensions,NPY_DOUBLE);
    if (!tmp_array)
    {
      cleanup();
      throw PythonException();
    }

  // PyObject argument is not an int ... try to build a contiguous PyArrayObject from it
  }
  else
  {
    if (!tmp_array)
      tmp_array = (PyArrayObject *)
	PyArray_ContiguousFromObject(pyObject,NPY_DOUBLE,0,0);

    // If this fails, clean up and throw a PythonException
    if (!tmp_array)
    {
      cleanup();
      throw PythonException();
    }
    // If the contiguous PyArrayObject built successfully, make sure it has the correct
    // number of dimensions
    else
    {
      int nd = PyArray_NDIM(tmp_array);
      npy_intp dimensions[ ] = { 1, blockMap.NumMyPoints() };  // Default dimensions
      bool reallocate = false;
      npy_intp arraySize = PyArray_MultiplyList(PyArray_DIMS(tmp_array),nd);

      if (nd < 2)
      {
	reallocate = true;
      }
      else
      {
	arraySize /= PyArray_DIMS(tmp_array)[0];
	if (arraySize != dimensions[1])
	{
	  dimensions[0] = PyArray_DIMS(tmp_array)[0];
	  reallocate = true;
	}
      }

      // Reallocate the tmp_array if necessary
      if (reallocate)
      {
	PyArrayObject * myArray = (PyArrayObject *)
	  PyArray_SimpleNew(2,dimensions,NPY_DOUBLE);
	if (!myArray)
	{
	  cleanup();
	  throw PythonException();
	}
	double        * myData  = (double *) PyArray_DATA(myArray);
	double        * tmpData = (double *) PyArray_DATA(tmp_array);
	for (int i=0; i<dimensions[0]; i++)
	{
	  for (int j=0; j<arraySize && j<dimensions[1]; j++)
	  {
	    myData[i*dimensions[1]+j] = tmpData[i*arraySize+j];
	  }
	}
	Py_XDECREF(tmp_array);
	tmp_array = myArray;
      }
    }
  }
  return (double *) (PyArray_DATA(tmp_array));
}

// =============================================================================
Epetra_Map & Epetra_NumPyMultiVector::getMap(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  if (!tmp_map)
  {
    const int totalLength = PyArray_Size((PyObject *)tmp_array);
    const int numVectors  = PyArray_DIMS(tmp_array)[0];
    tmp_map = new Epetra_Map(totalLength/numVectors,0,defaultComm);
  }
  return *tmp_map;
}

// =============================================================================
int Epetra_NumPyMultiVector::getNumVectors(PyObject * pyObject)
{
  if (!tmp_array) getArray(pyObject);
  return PyArray_DIMS(tmp_array)[0];
}

// =============================================================================
int Epetra_NumPyMultiVector::getNumVectors(const Epetra_BlockMap & blockMap,
					   PyObject * pyObject)
{
  if (!tmp_array) getArray(blockMap,pyObject);
  return PyArray_DIMS(tmp_array)[0];
}

// =============================================================================
int Epetra_NumPyMultiVector::getVectorSize(PyObject * pyObject)
{
  if (!tmp_map) getMap(pyObject);
  return tmp_map->NumMyPoints();
}

// =============================================================================
int * Epetra_NumPyMultiVector::getRange(PyObject * range,
					const Epetra_MultiVector & source)
{
  // Handle the default case (range == NULL), which is to return a
  // range of all the Epetra_MultiVector vectors
  if (range == NULL)
  {
    npy_intp dims[ ] = { (npy_intp) source.NumVectors() };
    tmp_range = (PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_INT);
    if (!tmp_range)
    {
      cleanup();
      throw PythonException();
    }
    int * data = (int *) PyArray_DATA(tmp_range);
    for (int i=0; i<dims[0]; i++) data[i] = i;
  }

  // Try to create a contiguous array of integers from the PyObject
  if (!tmp_range)
    tmp_range = (PyArrayObject *)
      PyArray_ContiguousFromObject(range,NPY_INT,1,1);

  // If this fails, clean up and throw a PythonException
  if (!tmp_range)
  {
    cleanup();
    throw PythonException();
  }

  // Obtain the length and return the array of integers
  return (int *) (PyArray_DATA(tmp_range));
}

// =============================================================================
int Epetra_NumPyMultiVector::getRangeLen(PyObject * range,
					 const Epetra_MultiVector & source)
{
  if (!tmp_range) getRange(range, source);
  return PyArray_Size((PyObject*)tmp_range);
}

// =============================================================================
void Epetra_NumPyMultiVector::cleanup()
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
  if (tmp_range)
  {
    Py_DECREF(tmp_range);
    tmp_range = NULL;
  }
}

// =============================================================================

// Constructors
// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap,
						 int numVectors, bool zeroOut):
  Epetra_MultiVector(blockMap, numVectors, zeroOut)
{
  // Create the array object
  npy_intp dims[ ] = { numVectors, blockMap.NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,
						      (void *)v[0]);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
  // Copy the Epetra_BlockMap
  map = new Epetra_BlockMap(blockMap);
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(const Epetra_MultiVector & source):
  Epetra_MultiVector(source)
{
  map = new Epetra_BlockMap(source.Map());
  npy_intp dims[ ] = { NumVectors(), map->NumMyPoints() };
  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,
						      (void *)v[0]);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }
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
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(Epetra_DataAccess CV,
						 const Epetra_NumPyMultiVector & source,
						 PyObject * range):
  Epetra_MultiVector(CV, (const Epetra_MultiVector &) source,
		     getRange(   range, (const Epetra_MultiVector) source),
		     getRangeLen(range, (const Epetra_MultiVector) source))
{
  // Store the Epetra_NumPyMultiVector's map
  map = new Epetra_BlockMap(source.Map());

  // Inintialize the local numpy array
  PyArrayObject * src_array = (PyArrayObject *) (source.ExtractView());
  int nd;
  // This shouldn't happen, but it does . . .
  if (NULL == src_array) nd = 2;
  else nd = PyArray_NDIM(src_array);
  npy_intp * dims = new npy_intp[nd];
  dims[0] = NumVectors();
  if (NULL == src_array) dims[1] = source.MyLength();
  else for (int i=1; i<nd; i++) dims[i] = PyArray_DIMS(src_array)[i];

  double **v = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_SimpleNewFromData(nd,dims,NPY_DOUBLE,
						      (void *)v[0]);
  delete [] dims;
  if (!array)
  {
    cleanup();
    throw PythonException();
  }

  // We're done with the tmp_range array
  Py_XDECREF(tmp_range);
  tmp_range = NULL;
}

// =============================================================================
Epetra_NumPyMultiVector::Epetra_NumPyMultiVector(Epetra_DataAccess CV,
						 const Epetra_MultiVector & source,
						 PyObject * range):
  Epetra_MultiVector(CV, source, getRange(range, source), getRangeLen(range, source))
{
  // Store the local map
  map = new Epetra_BlockMap(source.Map());

  // Wrap the Epetra_MultiVector
  npy_intp dims[ ] = { NumVectors(), MyLength() };
  double **v  = NULL;
  Epetra_MultiVector::ExtractView(&v);
  array = (PyArrayObject *) PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,
						      (void *)v[0]);
  if (!array)
  {
    cleanup();
    throw PythonException();
  }

  // We're done with the tmp_range array
  Py_XDECREF(tmp_range);
  tmp_range = NULL;
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
}

// =============================================================================
// Destructor
Epetra_NumPyMultiVector::~Epetra_NumPyMultiVector()
{
  Py_XDECREF(array);
  delete map;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::ExtractCopy() const
{
  return PyArray_NewCopy(array,NPY_ANYORDER);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::ExtractView() const
{
  Py_INCREF(array);
  return PyArray_Return(array);
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Dot(const Epetra_MultiVector & a) const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::Dot(a, result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Dot returned error code %d", status);
    goto fail;
  }
  po = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA(((PyArrayObject*)po));
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm1() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::Norm1(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Norm1 returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::Norm2() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::Norm2(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "Norm2 returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::NormInf() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;  

  int status = Epetra_MultiVector::NormInf(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "NormInf returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::
NormWeighted(const Epetra_MultiVector & weights) const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;  

  int status = Epetra_MultiVector::NormWeighted(weights,result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "NormWeighted returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MinValue() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::MinValue(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MinValue returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MaxValue() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::MaxValue(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MaxValue returned error code %d", status);
    goto fail;
  }
  po   = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

// =============================================================================
PyObject * Epetra_NumPyMultiVector::MeanValue() const
{
  int n = NumVectors();
  double * result = new double[n];
  npy_intp numVectors[ ] = {n};
  PyObject * po;
  double * data;

  int status = Epetra_MultiVector::MeanValue(result);

  if (status)
  {
    PyErr_Format(PyExc_RuntimeError, "MeanValue returned error code %d", status);
    goto fail;
  }
  po = PyArray_SimpleNew(1, numVectors, NPY_DOUBLE);
  data = (double*) PyArray_DATA((PyArrayObject*)po);
  for (int i=0; i<n; i++) data[i] = result[i];
  delete [] result;
  return PyArray_Return((PyArrayObject*)po);
 fail:
  delete [] result;
  return NULL;
}

}  // Namespace PyTrilinos
