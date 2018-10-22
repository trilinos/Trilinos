// -*- c++ -*-

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

%define %pliris_docstring
"
PyTrilinos.Pliris is the python interface to the Trilinos package
Pliris, an LU solver for dense matrices:

    http://trilinos.sandia.gov/packages/pliris

The purpose of Pliris is to provide an object-oriented interface to an
LU solver for dense matrices on parallel platforms. These matrices are
double precision real matrices distributed on a parallel machine.  The
python version of the Pliris package supports the following class:

    * Pliris  - Primary solver class.
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %pliris_docstring) Pliris

%{
// PyTrilinos configuration
#include "PyTrilinos_config.h"

// Epetra include files
#ifdef HAVE_EPETRA
// #undef HAVE_INTTYPES_H
// #undef HAVE_STDINT_H
#include "PyTrilinos_Epetra_Headers.hpp"

// NumPy include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Pliris include files
#include "PyTrilinos_Pliris_Headers.hpp"
#endif
%}

// Auto-documentation feature
%feature("autodoc", "1");

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "std_except.i"
%include "std_string.i"
using std::string;

// Include Pliris documentation
// %include "Pliris_dox.i"

// External Trilinos modules
#ifdef HAVE_EPETRA
%ignore Epetra_Version();
%import "Epetra.i"
#endif

///////////////////////////////////
// Pliris configuration support //
///////////////////////////////////
%include "Pliris_config.h"

////////////////////////////
// Pliris Version support //
////////////////////////////
%rename(Version) Pliris_Version;
%include "PlirisVersion.h"
%pythoncode
{
  __version__ = Version()
}

////////////////////
// Pliris support //
////////////////////
%include "Pliris.h"
