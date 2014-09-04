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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYTRILINOS_DOMI_UTIL_HPP
#define PYTRILINOS_DOMI_UTIL_HPP

// Include the PyTrilinos Distributed Array Protocol header
#include "PyTrilinos_DAP.hpp"

// Include Domi headers
#include "Domi_MDMap.hpp"
#include "Domi_MDVector.hpp"

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

// Convert a Domi::MDArrayRCP to a NumPy array.  It is assumes that
// the user will call this function with the NumPy typecode that
// corresponds to the template parameter T.
template< class T >
PyObject *
convertToNumPyArray(const Domi::MDArrayRCP< T > & mdArrayRcp,
                    int typecode);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDComm.
Teuchos::RCP< const Domi::MDComm >
convertToMDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDMap.
Teuchos::RCP< const Domi::MDMap<> >
convertToMDMap(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDVector.
template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar > >
convertToMDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                  const DistArrayProtocol & distarray);

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
convertToNumPyArray(const Domi::MDArrayRCP< T > & mdArrayRcp,
                    int typecode)
{
  // Get the number of dimensions and initialize the dimensions and
  // strides arrays
  int numDims = mdArrayRcp.numDims();
  Teuchos::Array< npy_intp > dims(numDims);
  Teuchos::Array< npy_intp > strides(numDims);

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

template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar > >
convertToMDVector(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                  const DistArrayProtocol & distarray)
{
  // Get the equivalent MDMap
  Teuchos::RCP< const Domi::MDMap<> > mdMap =
    convertToMDMap(teuchosComm, distarray);

  // Get the equivalent MDArrayRCP
  Domi::MDArrayRCP< Scalar > mdArrayRcp =
    convertToMDArrayRCP< Scalar >((PyArrayObject*) distarray.buffer());

#ifdef PYTRILINOS_DOMI_UTIL_VERBOSE
  std::cout << "mdArrayRcp = " << mdArrayRcp << std::endl;
#endif

  // Return the result
  return Teuchos::rcp(new Domi::MDVector< Scalar >(mdMap, mdArrayRcp));
}

////////////////////////////////////////////////////////////////////////

}

#endif
