// -*- c++ -*-

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

// Domi_MDArray.i is a SWIG interface file that provides SWIG
// directives to handle Domi MDArray types (including MDArrayView and
// MDArrayRCP).  These classes are not wrapped, but instead typemaps
// are defined so that the python user can provide NumPy arrays as
// input, or accept NumPy arrays as output, instead.  Currently, the
// following classes are handled:
//
//     Domi::MDArray< T >
//     Domi::MDArrayView< T >
%{
#include "Domi_MDArrayView.hpp"
using Domi::MDArrayView;
%}
#define REFCOUNTPTR_INLINE
%import  "Domi_MDArrayView.hpp"

////////////////////////////////////////////////////////////////////////
// The philosophy is that wherever Domi MDArray classes are used in
// C++, NumPy arrays will be used in python.  Thus we need the NumPy
// SWIG directives.
%include "numpy.i"

////////////////////////////////////////////////////////////////////////
// Define a macro that takes a C++ data type (TYPE) and a
// corresponding NumPy typecode (TYPECODE) and define all of the
// typemaps needed to handle that TYPE array.
%define %domi_mdarray_typemaps(TYPE, TYPECODE)

// If an MDArrayView argument has a const TYPE, then we know that the
// argument is input only.  Therefore we allow any type of sequence to
// be converted to a PyArrayObject and then extract the resulting data
// pointer to construct the MDArrayView.  If the conversion creates a
// new PyArrayObject, then we have to be sure to decrement its
// reference count once the MDArrayView has been used.
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Domi::MDArrayView< const TYPE >)
{
  $1 = is_array($input) || PyList_Check($input) || PyTuple_Check($input);
}

%typemap(in) Domi::MDArrayView< const TYPE >
(int is_new = 0,
 PyArrayObject * npArray = NULL)
{
  npArray = obj_to_array_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  Domi::size_type npSize = 1;
  int numDims = array_numdims(npArray);
  Teuchos::Array< Domi::dim_type  > dims(numDims);
  Teuchos::Array< Domi::size_type > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis]  = array_size(  npArray, axis);
    strides[axis]  = array_stride(npArray, axis);
    npSize        += dims[axis] * strides[axis];
    strides[axis] *= sizeof(TYPE);
  }
  Teuchos::ArrayView< TYPE > data((TYPE*) array_data(npArray), npSize);
  Domi::Layout layout = PyArray_IS_C_CONTIGUOUS(npArray) ?
                        Domi::C_ORDER : Domi::FORTRAN_ORDER;
  $1 = Domi::MDArrayView< TYPE >(data, dims, strides, layout);
}

%typemap(freearg) Domi::MDArrayView< const TYPE >
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Domi::MDArrayView< TYPE > const &)
{
  $1 = is_array($input) || PyList_Check($input) || PyTuple_Check($input);
}

%typemap(in) Domi::MDArrayView< TYPE > const &
(int is_new = 0,
 PyArrayObject * npArray = NULL,
 Domi::MDArrayView< TYPE > temp)
{
  npArray = obj_to_array_allow_conversion($input, TYPECODE, &is_new);
  if (!npArray) SWIG_fail;
  Domi::size_type npSize = 1;
  int numDims = array_numdims(npArray);
  Teuchos::Array< Domi::dim_type  > dims(numDims);
  Teuchos::Array< Domi::size_type > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis]  = array_size(  npArray, axis);
    strides[axis]  = array_stride(npArray, axis);
    npSize        += dims[axis] * strides[axis];
    strides[axis] *= sizeof(TYPE);
  }
  Teuchos::ArrayView< TYPE > data((TYPE*) array_data(npArray), npSize);
  Domi::Layout layout = PyArray_IS_C_CONTIGUOUS(npArray) ?
                        Domi::C_ORDER : Domi::FORTRAN_ORDER;
  temp = Domi::MDArrayView< TYPE >(data, dims, strides, layout);
  $1 = &temp;
}

%typemap(freearg) Domi::MDArrayView< TYPE > const &
{
  if (is_new$argnum) Py_DECREF(npArray$argnum);
}

// If an MDArrayView argument has a non-const TYPE, then the default
// behavior is to assume that the array is input/output.  Therefore
// the input python argument must be a NumPy array.
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (Domi::MDArrayView< TYPE >)
{
  $1 = is_array($input);
}

%typemap(in) Domi::MDArrayView< TYPE >
{
  PyArrayObject * npArray = obj_to_array_no_conversion($input, TYPECODE);
  if (!npArray) SWIG_fail;
  Domi::size_type npSize = 1;
  int numDims = array_numdims(npArray);
  Teuchos::Array< Domi::dim_type  > dims(numDims);
  Teuchos::Array< Domi::size_type > strides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dims[   axis]  = array_size(  npArray, axis);
    strides[axis]  = array_stride(npArray, axis);
    npSize        += dims[axis] * strides[axis];
    strides[axis] *= sizeof(TYPE);
  }
  Teuchos::ArrayView< TYPE > data((TYPE*) array_data(npArray), npSize);
  Domi::Layout layout = PyArray_IS_C_CONTIGUOUS(npArray) ?
                        Domi::C_ORDER : Domi::FORTRAN_ORDER;
  $1 = Domi::MDArrayView< TYPE >(data, dims, strides, layout);
}

// If an MDArray or MDArrayView is output, with either a const or
// non-const TYPE, convert the underlying data to a NumPy array of
// correct type.
%typemap(out) Domi::MDArray< TYPE >
{
  int numDims = $1.numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1.dimension(axis);
    npStrides[axis] = $1.strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1.layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY : NPY_ARRAY_FARRAY;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1.getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%typemap(out) Domi::MDArrayView< TYPE >
{
  int numDims = $1.numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1.dimension(axis);
    npStrides[axis] = $1.strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1.layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY : NPY_ARRAY_FARRAY;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1.getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%typemap(out) Domi::MDArray< TYPE > const &
{
  int numDims = $1->numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1->dimension(axis);
    npStrides[axis] = $1->strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1->layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY_RO : NPY_ARRAY_FARRAY_RO;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1->getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%typemap(out) Domi::MDArrayView< TYPE > const &
{
  int numDims = $1->numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1->dimension(axis);
    npStrides[axis] = $1->strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1->layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY_RO : NPY_ARRAY_FARRAY_RO;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1->getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%typemap(out) Domi::MDArray< const TYPE >
{
  int numDims = $1.numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1.dimension(axis);
    npStrides[axis] = $1.strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1.layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY_RO : NPY_ARRAY_FARRAY_RO;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1.getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%typemap(out) Domi::MDArrayView< const TYPE >
{
  int numDims = $1.numDims();
  Teuchos::Array< npy_intp > npDims(   numDims);
  Teuchos::Array< npy_intp > npStrides(numDims);
  for (int axis = 0; axis < numDims; ++axis)
  {
    npDims[   axis] = $1.dimension(axis);
    npStrides[axis] = $1.strides()[axis] * sizeof(TYPE);
  }
  int flags = ($1.layout() == Domi::C_ORDER) ?
              NPY_ARRAY_CARRAY_RO : NPY_ARRAY_FARRAY_RO;
  $result = PyArray_New(&PyArray_Type,
                        numDims,
                        npDims.getRawPtr(),
                        TYPECODE,
                        npStrides.getRawPtr(),
                        (void*) $1.getRawPtr(),
                        -1, //sizeof(TYPE),
                        flags,
                        NULL);
  if (!$result) SWIG_fail;
}

%enddef

////////////////////////////////////////////////////////////////////////
// Call the %domi_mdarray_typemaps() macro for specific data types
// that are supported by NumPy
%domi_mdarray_typemaps(signed char       , NPY_BYTE     )
%domi_mdarray_typemaps(unsigned char     , NPY_UBYTE    )
%domi_mdarray_typemaps(short             , NPY_SHORT    )
%domi_mdarray_typemaps(unsigned short    , NPY_USHORT   )
%domi_mdarray_typemaps(int               , NPY_INT      )
%domi_mdarray_typemaps(unsigned int      , NPY_UINT     )
%domi_mdarray_typemaps(long              , NPY_LONG     )
%domi_mdarray_typemaps(unsigned long     , NPY_ULONG    )
%domi_mdarray_typemaps(long long         , NPY_LONGLONG )
%domi_mdarray_typemaps(unsigned long long, NPY_ULONGLONG)
%domi_mdarray_typemaps(float             , NPY_FLOAT    )
%domi_mdarray_typemaps(double            , NPY_DOUBLE   )
