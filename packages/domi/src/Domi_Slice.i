// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2014) Sandia Corporation
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

%{
// Domi includes
#include "Domi_Slice.hpp"

// PyTrilinos includes
#include "PyTrilinos_Domi_Util.hpp"
%}

// Type checking for Domi::Slice.  We do not wrap the Slice class, and
// instead only accept Python slices as input
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) Domi::Slice
{
  $1 = PySlice_Check($input);
}
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) const Domi::Slice
{
  $1 = PySlice_Check($input);
}
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) Domi::Slice &
{
  $1 = PySlice_Check($input);
}
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) const Domi::Slice &
{
  $1 = PySlice_Check($input);
}

// Input typemaps for Domi::Slice.  We only accept Python slices as
// input, and so we convert directly to Domi::Slice.  The C/API
// requires that we provide the length of the container, which is not
// available to us in this general setting, so we provide the
// Domi::Slice::Default value, which will give us a Domi::Slice with
// Default values if the Python slice has default values.
%typemap(in) Domi::Slice
{
  $1 = PyTrilinos::convertToDomiSlice((PySliceObject*) $input,
                                      (Py_ssize_t) Domi::Slice::Default);
}
%typemap(in) const Domi::Slice
{
  $1 = PyTrilinos::convertToDomiSlice((PySliceObject*) $input,
                                      (Py_ssize_t) Domi::Slice::Default);
}
%typemap(in) Domi::Slice &
(Domi::Slice tempSlice)
{
  tempSlice =
    PyTrilinos::convertToDomiSlice((PySliceObject*) $input,
                                   (Py_ssize_t) Domi::Slice::Default);
  $1 = &tempSlice;
}
%typemap(in) const Domi::Slice &
(Domi::Slice tempSlice)
{
  tempSlice =
    PyTrilinos::convertToDomiSlice((PySliceObject*) $input,
                                   (Py_ssize_t) Domi::Slice::Default);
  $1 = &tempSlice;
}

// When a Domi::Slice object is returned from a C++ method, convert it
// to a Python slice.  Default Domi::Slice values will be converted to
// default Python slice values
%typemap(out) Domi::Slice
{
  $result = PyTrilinos::convertToPySlice($1);
}
