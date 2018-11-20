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

%define %ifpack_docstring
"
PyTrilinos.IFPACK is the python interface to the Trilinos
preconditioner package IFPACK:

    http://trilinos.sandia.gov/packages/ifpack

The purpose of IFPACK is to provide incomplete foctorization
preconditioners to Trilinos.  Note that the C++ version of IFPACK uses
the prefix 'Ifpack_' which has been stripped from the python version.

IFPACK provides the following user-level classes:

    * Factory         - A factory for creating IFPACK preconditioners
    * Preconditioner  - Pure virtual base class for defining interface
    * IC              - Incomplete Cholesky preconditioner
    * ICT             - Incomplete Cholesky preconditioner w/threshold
    * ILU             - Incomplete lower/upper preconditioner
    * ILUT            - Incomplete lower/upper preconditioner w/threshold
    * PointRelaxation - Point relaxation predonditioner
    * Amesos          - Use Amesos factorizations as preconditioners

and functions:

    * AnalyzeMatrix          - Analyze the basic properties of a matrix
    * AnalyzeMatrixElements  - Analyze the distribution of values of a matrix
    * AnalyzeVectorElements  - Analyze the distribution of values of a vector
    * PrintSparsity          - Create PS file with sparsity pattern of matrix

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exIFPACK.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %ifpack_docstring) IFPACK

%{
// System include files
#include <iostream>
#include <sstream>
#include <vector>

// Configuration include files
#include "PyTrilinos_config.h"

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"

// NumPy include
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#endif

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#include "PyTrilinos_Epetra_Headers.hpp"

// IFPACK include files
#include "PyTrilinos_IFPACK_Headers.hpp"
%}

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include IFPACK documentation
%include "IFPACK_dox.i"

// External Trilinos modules
%import "Teuchos.i"
#ifdef HAVE_PYTRILINOS_EPETRA
%ignore Epetra_Version();
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
    SWIG_exception(SWIG_UnknownError, "Unkown C++ exception");
  }
}

// General ignore directives
%ignore operator<<;

////////////////////////
// I/O Stream support //
////////////////////////
// #pragma SWIG nowarn=401
//%include "std_iostream.i"

//////////////////////////////////
// IFPACK configuration support //
//////////////////////////////////
#undef  PACKAGE_BUGREPORT
%ignore PACKAGE_BUGREPORT;
#undef  PACKAGE_NAME
%ignore PACKAGE_NAME;
#undef  PACKAGE_STRING
%ignore PACKAGE_STRING;
#undef  PACKAGE_TARNAME
%ignore PACKAGE_TARNAME;
#undef  PACKAGE_VERSION
%ignore PACKAGE_VERSION;
%include "Ifpack_config.h"
%include "Ifpack_ConfigDefs.h"

////////////////////////////
// IFPACK_Version support //
////////////////////////////
%include "Ifpack_Version.h"
%pythoncode
%{
Version = Ifpack_Version
__version__ = Version().split()[3]
%}

//////////////////////////
// IFPACK_Utils support //
//////////////////////////
%rename(AnalyzeMatrix        ) Ifpack_Analyze;
%rename(AnalyzeMatrixElements) Ifpack_AnalyzeMatrixElements;
%rename(AnalyzeVectorElements) Ifpack_AnalyzeVectorElements;
%rename(PrintSparsity        ) Ifpack_PrintSparsity;
%include "Ifpack_Utils.h"

///////////////////////////////////
// IFPACK_Preconditioner support //
///////////////////////////////////
%teuchos_rcp(Ifpack_Preconditioner)
%ignore operator<<(ostream &, const Ifpack_Preconditioner &);
%ignore Ifpack_Preconditioner::Condest() const;
%rename(Preconditioner) Ifpack_Preconditioner;
%include "Ifpack_Preconditioner.h"
%extend Ifpack_Preconditioner
{
  std::string __str__()
  {
    std::stringstream os;
    os << *self;
    return os.str();
  }
  //void __del__()
  //{
  //  delete self;
  //}
}

///////////////////////
// IFPACK_IC support //
///////////////////////
%teuchos_rcp(Ifpack_IC)
%ignore Ifpack_IC::Condest() const;
%rename(IC) Ifpack_IC;
%include "Ifpack_IC.h"

////////////////////////
// IFPACK_ICT support //
////////////////////////
%teuchos_rcp(Ifpack_ICT)
%ignore Ifpack_ICT::Condest() const;
%rename(ICT) Ifpack_ICT;
%include "Ifpack_ICT.h"

////////////////////////
// IFPACK_ILU support //
////////////////////////
%teuchos_rcp(Ifpack_ILU)
%ignore Ifpack_ILU::Condest() const;
%rename(ILU) Ifpack_ILU;
%include "Ifpack_ILU.h"

/////////////////////////
// IFPACK_ILUT support //
/////////////////////////
%teuchos_rcp(Ifpack_ILUT)
%ignore Ifpack_ILUT::Condest() const;
%rename(ILUT) Ifpack_ILUT;
%include "Ifpack_ILUT.h"

////////////////////////////////////
// IFPACK_PointRelaxation support //
////////////////////////////////////
%teuchos_rcp(Ifpack_PointRelaxation)
%ignore Ifpack_PointRelaxation::Condest() const;
%rename(PointRelaxation) Ifpack_PointRelaxation;
%include "Ifpack_PointRelaxation.h"

///////////////////////////
// IFPACK_Amesos support //
///////////////////////////
%teuchos_rcp(Ifpack_Amesos)
%ignore Ifpack_Amesos::Condest() const;
%rename(Amesos) Ifpack_Amesos;
%include "Ifpack_Amesos.h"

////////////////////////////
// IFPACK factory support //
////////////////////////////
%ignore operator>>(std::istream &, Ifpack::EPrecType &);
%newobject Ifpack::Create;
%rename(Factory) Ifpack;
%include "Ifpack.h"

// Turn off the exception handling
%exception;
