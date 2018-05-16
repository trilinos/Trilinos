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

%define %triutils_docstring
"
PyTrilinos.TriUtils is the python interface to the Trilinos utilities
package TriUtils:

    http://trilinos.sandia.gov/packages/triutils

The purpose of TriUtils is to provide some utilities typically needed
when testing Trilinos software.

TriUtils provides the following user-level classes:

    * CrsMatrixGallery  - Provide example CrsMatrix objects
    * VbrMatrixGallery  - Provide example VbrMatrix objects

and function:

    * ReadHB            - Obtain a problem from an HB file

For an examples of usage, please consult the following script in the
example subdirectory of the PyTrilinos package:

    * exIFPACK.py
"
%enddef

%module(package = "PyTrilinos",
	autodoc = "1",
	docstring = %triutils_docstring) TriUtils

%{
// System include files
#include <iostream>
#include <sstream>
#include <vector>

// Configuration include files
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Trilinos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"

// NumPy include
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#endif

// Trilinos utility include files
#include "PyTrilinos_TriUtils_Headers.hpp"
%}

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Auto-documentation feature
%feature("autodoc", "1");

// Include the TriUtils documentation
%include "TriUtils_dox.i"    // Doxygen-generated documentation

// Standard exception handling
%include "exception.i"

// General ignore directives
#pragma SWIG nowarn=503
%ignore *::operator<< ;

// Epetra interface include files
#ifdef HAVE_EPETRA
%import "Epetra.i"
#endif

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}


%exception
{
  try
  {
    $action
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

///////////////////////////////////
// Trilinos_Util_Version support //
///////////////////////////////////
%rename (TriUtils_Version) Triutils_Version;
%include "Trilinos_Util_Version.h"
%pythoncode
%{
__version__ = TriUtils_Version().split()[3]
%}

///////////////////////////
// Trilinos_Util support //
///////////////////////////
#ifdef HAVE_EPETRA
%ignore Trilinos_Util_read_hb;
%ignore Trilinos_Util_read_coo;
%rename (ReadHB) Trilinos_Util_ReadHb2Epetra;
%rename (ReadHB64) Trilinos_Util_ReadHb2Epetra64;
%rename (ReadHpc) Trilinos_Util_ReadHpc2Epetra;
%rename (ReadHpc64) Trilinos_Util_ReadHpc2Epetra;
%rename (ReadHBVbr) Trilinos_Util_ReadHb2EpetraVbr;
%rename (ReadHBVbr64) Trilinos_Util_ReadHb2EpetraVbr64;
%ignore Trilinos_Util_distrib_msr_matrix;
%ignore Trilinos_Util_distrib_vbr_matrix;
%ignore Trilinos_Util_create_vbr;
%ignore Trilinos_Util_smsrres;
%ignore Trilinos_Util_scscres;
%ignore Trilinos_Util_scscmv;
%ignore Trilinos_Util_svbrres;
%ignore Trilinos_Util_msr2vbr;
%ignore Trilinos_Util_find_block_col;
%ignore Trilinos_Util_find_block_in_row;
%ignore Trilinos_Util_add_new_ele;
%ignore Trilinos_Util_find_closest_not_larger;
%ignore Trilinos_Util_convert_values_to_ptrs;
%ignore Trilinos_Util_csrcsc;
%ignore Trilinos_Util_csrmsr;
%ignore Trilinos_Util_ssrcsr;
%ignore Trilinos_Util_coocsr;
%ignore SPBLASMAT_STRUCT;
%ignore SPBLASMAT;
%ignore Trilinos_Util_duscr_vbr;
%ignore Trilinos_Util_dusmm;
%ignore Trilinos_Util_dusds_vbr;
%rename (GenerateCrsProblem) Trilinos_Util_GenerateCrsProblem;
%rename (GenerateCrsProblem64) Trilinos_Util_GenerateCrsProblem64;
%rename (GenerateVbrProblem) Trilinos_Util_GenerateVbrProblem;
%rename (ReadTriples) Trilinos_Util_ReadTriples2Epetra;
%rename (ReadTriples64) Trilinos_Util_ReadTriples2Epetra64;
%ignore Trilinos_Util_write_vec;
%ignore Trilinos_Util_read_vec;
%include "Trilinos_Util.h"
#endif

////////////////////////////////////////////
// Trilinos_Util_CrsMatrixGallery support //
////////////////////////////////////////////
%ignore
Trilinos_Util::CrsMatrixGallery::operator<<(ostream&,
					    const Trilinos_Util::CrsMatrixGallery&);
%include "Trilinos_Util_CrsMatrixGallery.h"

// Turn off the exception handling
%exception;
