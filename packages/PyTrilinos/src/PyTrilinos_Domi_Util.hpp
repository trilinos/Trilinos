// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
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

#ifndef PYTRILINOS_DOMI_UTIL_HPP
#define PYTRILINOS_DOMI_UTIL_HPP

// Include PyTrilinos utilities
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_NumPy_Util.hpp"
#include "PyTrilinos_DAP.hpp"

// Include Domi headers
#include "Domi_MDMap.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Domi_MDVector.hpp"
#pragma GCC diagnostic pop

// Verbosity
// #define PYTRILINOS_DOMI_UTIL_VERBOSE

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

// Convert a Domi::Slice to a Python slice object.  Any Domi::Slice
// values equal to Domi::Slice::Default will be converted to None
// values in the Python slice.
PyObject *
convertToPySlice(const Domi::Slice & domiSlice);

////////////////////////////////////////////////////////////////////////

// Convert a Python slice object to a Domi::Slice.  The Python C/API
// for slices requires the user to provide a length, so this converter
// does as well.  If the length is not known, then use
// Domi::Slice::Default.
Domi::Slice
convertToDomiSlice(PySliceObject * pySlice,
                   Py_ssize_t length);

////////////////////////////////////////////////////////////////////////

// Convert a NumPy array object to an MDArrayRCP.  It is assumed that
// the user has already checked the NumPy array's data type and is
// calling this function with the proper template parameter T.
template< class T >
Domi::MDArrayRCP< T >
convertToMDArrayRCP(PyArrayObject * pyArray);

////////////////////////////////////////////////////////////////////////

// Convert a Domi::MDArrayRCP to a NumPy array.
template< class T >
PyObject *
convertToNumPyArray(const Domi::MDArrayRCP< T > & mdArrayRcp);

////////////////////////////////////////////////////////////////////////

// Convert a Domi::MDArrayView to a NumPy array.
template< class T >
PyObject *
convertToNumPyArray(const Domi::MDArrayView< T > & mdArrayView);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDComm.
Teuchos::RCP< const Domi::MDComm >
convertToMDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDMap.
Teuchos::RCP< const Domi::MDMap >
convertToMDMap(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

PyObject * convertToDimData(const Teuchos::RCP< const Domi::MDMap > & mdMap);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDVector.
template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar > >
convertToMDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                  const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

template< class Scalar >
PyObject *
convertToDistArray(Domi::MDVector< Scalar > & mdVector);

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Implementation of templated functions
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////

template< class T >
Domi::MDArrayRCP< T >
convertToMDArrayRCP(PyArrayObject * pyArray)
{
  // Get the number of dimensions and initialize the dimensions and
  // strides arrays
  int numDims = PyArray_NDIM(pyArray);
  Teuchos::Array< Domi::dim_type  > dims(   numDims);
  Teuchos::Array< Domi::size_type > strides(numDims);

  // Set the dimensions and strides
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis] = (Domi::dim_type ) PyArray_DIM(   pyArray, axis);
    strides[axis] = (Domi::size_type) PyArray_STRIDE(pyArray, axis);
  }

  // Get the data pointer and layout
  T * data = (T*) PyArray_DATA(pyArray);
  Domi::Layout layout = PyArray_IS_C_CONTIGUOUS(pyArray) ? Domi::C_ORDER :
    Domi::FORTRAN_ORDER;

  // Return the result
  return Domi::MDArrayRCP< T >(dims, strides, data, layout);
}

////////////////////////////////////////////////////////////////////////

template< class T >
PyObject *
convertToNumPyArray(const Domi::MDArrayRCP< T > & mdArrayRcp)
{
  // Get the number of dimensions and initialize the dimensions and
  // strides arrays
  int numDims = mdArrayRcp.numDims();
  Teuchos::Array< npy_intp > dims(numDims);
  Teuchos::Array< npy_intp > strides(numDims);
  int typecode = NumPy_TypeCode< T >();

  // Set the dimensions and strides
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis] = mdArrayRcp.dimension(axis);
    strides[axis] = mdArrayRcp.strides()[axis];
  }
  
  // Get the data pointer and flags, based on data layout
  void * data = (void*) mdArrayRcp.getRawPtr();
  int flags = (mdArrayRcp.layout() == Domi::C_ORDER) ? NPY_ARRAY_CARRAY :
    NPY_ARRAY_FARRAY;

  // Return the result
  return PyArray_New(&PyArray_Type,
                     numDims,
                     dims.getRawPtr(),
                     typecode,
                     strides.getRawPtr(),
                     data,
                     0,
                     flags,
                     NULL);
}

////////////////////////////////////////////////////////////////////////

template< class T >
PyObject *
convertToNumPyArray(const Domi::MDArrayView< T > & mdArrayView)
{
  // Get the number of dimensions and initialize the dimensions and
  // strides arrays
  int numDims = mdArrayView.numDims();
  Teuchos::Array< npy_intp > dims(numDims);
  Teuchos::Array< npy_intp > strides(numDims);
  int typecode = NumPy_TypeCode< T >();

  // Set the dimensions and strides
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis] = mdArrayView.dimension(axis);
    strides[axis] = mdArrayView.strides()[axis];
  }
  
  // Get the data pointer and flags, based on data layout
  void * data = (void*) mdArrayView.getRawPtr();
  int flags = (mdArrayView.layout() == Domi::C_ORDER) ? NPY_ARRAY_CARRAY :
    NPY_ARRAY_FARRAY;

  // Return the result
  return PyArray_New(&PyArray_Type,
                     numDims,
                     dims.getRawPtr(),
                     typecode,
                     strides.getRawPtr(),
                     data,
                     0,
                     flags,
                     NULL);
}

////////////////////////////////////////////////////////////////////////

template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar > >
convertToMDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                  const DistArrayProtocol & distarray)
{
  // Get the equivalent MDMap
  Teuchos::RCP< const Domi::MDMap > mdMap =
    convertToMDMap(teuchosComm, distarray);

  // Get the equivalent MDArrayRCP
  Domi::MDArrayRCP< Scalar > mdArrayRcp =
    convertToMDArrayRCP< Scalar >((PyArrayObject*) distarray.buffer());

#ifdef PYTRILINOS_DOMI_UTIL_VERBOSE
  std::cout << "mdArrayRcp = " << mdArrayRcp << std::endl;
#endif

  // Return the result
  try
  {
    return Teuchos::rcp(new Domi::MDVector< Scalar >(mdMap, mdArrayRcp));
  }
  catch (Domi::InvalidArgument & e)
  {
    PyErr_SetString(PyExc_ValueError, e.what());
    throw PythonException();
  }
}

////////////////////////////////////////////////////////////////////////

template< class Scalar >
PyObject *
convertToDistArray(Domi::MDVector< Scalar > & mdVector)
{
  PyObject * distArrayProtocol;
  PyObject * buffer;
  PyObject * dimData;

  distArrayProtocol = PyDict_New();
  if (!distArrayProtocol) goto fail;
  if (PyDict_SetItemString(distArrayProtocol, "__version__",
                           Py_BuildValue("s","0.10.0")) == -1) goto fail;

  buffer = convertToNumPyArray(mdVector.getDataNonConst());
  if (!buffer) goto fail;
  if (PyDict_SetItemString(distArrayProtocol, "buffer", buffer) == -1)
    goto fail;

  dimData = convertToDimData(mdVector.getMDMap());
  if (!dimData) goto fail;
  if (PyDict_SetItemString(distArrayProtocol, "dim_data", dimData) == -1)
    goto fail;

  Py_DECREF(buffer);
  Py_DECREF(dimData);
  return distArrayProtocol;

  fail:
  if (distArrayProtocol) PyDict_Clear(distArrayProtocol);
  Py_XDECREF(distArrayProtocol);
  Py_XDECREF(buffer);
  Py_XDECREF(dimData);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

}

#endif

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

