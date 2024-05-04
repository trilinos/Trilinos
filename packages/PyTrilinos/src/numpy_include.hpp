// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
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

#ifndef NUMPY_INCLUDE_HPP
#define NUMPY_INCLUDE_HPP

// This include file takes care of five of the six things necessary
// when including the numpy header file arrayobject.h.  First, the
// Python.h header file is included.  Second, the
// PY_ARRAY_UNIQUE_SYMBOL is defined, which will allow PyTrilinos to
// work with other extension modules that are compiled against NumPy.
// Third, the NPY_NO_DEPRECATED macro is set to NPY_1_7_API_VERSION to
// ensure that no deprecated NumPy code is used.  Fourth, the
// numpy/arrayobject.h header file is included.  Fifth and finally,
// the NPY_API_VERSION macro from arrayobject.h is checked, and if it
// is old enough, macros are defined so that PyTrilinos will compile
// with older versions of NumPy.

// The user is responsible for defining the macro NO_IMPORT_ARRAY in
// those source files that do not call the numpy routine
// import_array().

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL PyTrilinos_NumPy
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if NPY_API_VERSION < 0x00000007
#define NPY_ANYORDER         PyArray_ANYORDER
#define NPY_DOUBLE           PyArray_DOUBLE
#define NPY_INT              PyArray_INT
#define NPY_ARRAY_FARRAY     NPY_FARRAY
#define NPY_ARRAY_DEFAULT    NPY_DEFAULT
#define NPY_ARRAY_NOTSWAPPED NPY_NOTSWAPPED
#endif

#endif // NUMPY_INCLUDE_HPP

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

