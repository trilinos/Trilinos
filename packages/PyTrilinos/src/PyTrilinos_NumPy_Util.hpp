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

#ifndef PYTRILINOS_NUMPY_UTIL_HPP
#define PYTRILINOS_NUMPY_UTIL_HPP

// System include
#include <complex>

// NumPy include
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarraytypes.h>

// Note: I first attempted to '#include <numpy/arrayobject.h>', but
// simply including it, but not using it, led to a segmentation fault
// in the resulting wrapper file.  What I need are the NPY_TYPES
// enumeration definition, and including just ndarraytypes.h seems to
// do the trick.

namespace PyTrilinos
{

/////////////////////////////////////////////////
// Prototype for the NumPy type code converter //
/////////////////////////////////////////////////
template< typename TYPE >
int NumPy_TypeCode();

/////////////////////
// Specializations //
/////////////////////

template<>
int NumPy_TypeCode< bool >()
{
  return NPY_BOOL;
}

template<>
int NumPy_TypeCode< char >()
{
  return NPY_BYTE;
}

template<>
int NumPy_TypeCode< unsigned char >()
{
  return NPY_UBYTE;
}

template<>
int NumPy_TypeCode< short >()
{
  return NPY_SHORT;
}

template<>
int NumPy_TypeCode< unsigned short >()
{
  return NPY_USHORT;
}

template<>
int NumPy_TypeCode< int >()
{
  return NPY_INT;
}

template<>
int NumPy_TypeCode< unsigned int >()
{
  return NPY_UINT;
}

template<>
int NumPy_TypeCode< long >()
{
  return NPY_LONG;
}

template<>
int NumPy_TypeCode< unsigned long >()
{
  return NPY_ULONG;
}

template<>
int NumPy_TypeCode< long long >()
{
  return NPY_LONGLONG;
}

template<>
int NumPy_TypeCode< unsigned long long >()
{
  return NPY_ULONGLONG;
}

template<>
int NumPy_TypeCode< float >()
{
  return NPY_FLOAT;
}

template<>
int NumPy_TypeCode< double >()
{
  return NPY_DOUBLE;
}

template<>
int NumPy_TypeCode< long double >()
{
  return NPY_LONGDOUBLE;
}

template<>
int NumPy_TypeCode< std::complex< float > >()
{
  return NPY_CFLOAT;
}

template<>
int NumPy_TypeCode< std::complex< double > >()
{
  return NPY_CDOUBLE;
}

template<>
int NumPy_TypeCode< std::complex< long double > >()
{
  return NPY_CLONGDOUBLE;
}

}

#endif
