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

%{
// Domi include files
#include "Domi_Slice.hpp"

// PyTrilinos include files
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
