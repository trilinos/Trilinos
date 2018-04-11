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

%define %loca_epetra_interface_docstring
"
PyTrilinos.LOCA.Epetra.Interface is the python interface to namespace
Epetra::Interface of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra.Interface is to provide a concrete Epetra
implementation of LOCA interfaces.  The python version of
LOCA.Epetra.Interface supports the following classes:

    * Required                 - Provides a set of interfaces for users to
                                 provide information about the nonlinear
                                 problem to LOCA
    * MassMatrix               - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for the coefficients of
                                 time dependent terms
    * TimeDependent            - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for computing the shifted
                                 matrix
    * TimeDependentMatrixFree  - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for applying the shifted
                                 matrix in a matrix-free setting
"
%enddef

%module(package      = "PyTrilinos.LOCA.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_interface_docstring) Interface

%{
// NumPy include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// PyTrilinos include files
#include "PyTrilinos_PythonException.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#include "PyTrilinos_Epetra_Headers.hpp"

// NOX include files
#include "PyTrilinos_NOX_Epetra_Headers.hpp"

// LOCA::Epetra::Interface include files
#include "PyTrilinos_LOCA_Epetra_Headers.hpp"
%}

// General ignore directives
%ignore *::operator=;

// STL support
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Director exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PyTrilinos::PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  catch(int errCode)
  {
    PyErr_Format(PyExc_EpetraError, "Error code = %d\nSee stderr for details", errCode);
    SWIG_fail;
  }
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Learn about LOCA::Abstract::Iterator::StepStatus enumeration
%import "LOCA_Abstract_Iterator.H"

// Teuchos::RCPs typemaps
%teuchos_rcp(LOCA::Epetra::Interface::Required)
%teuchos_rcp(LOCA::Epetra::Interface::MassMatrix)
%teuchos_rcp(LOCA::Epetra::Interface::TimeDependent)
%teuchos_rcp(LOCA::Epetra::Interface::TimeDependentMatrixFree)

///////////////////////
// NOX_Utils support //
///////////////////////
%import "NOX_Utils.i"


%feature("autodoc",
"computeF(self, Epetra.Vector x, Epetra.Vector F, FillType flag) -> bool

  Virtual method in C++ that is intended to be overridden by user.
  This method defines the nonlinear function to be solved.  Arguments
  x and F will be provided as numpy-hybrid Epetra.Vector objects.
  Return True if the computation is successful.

  It is strongly recommended that computeF() not raise any exceptions,
  accidental or otherwise.  This can be prevented by wrapping your
  algorithm in a try block:

    try:
      # Your code here...
    except Exception, e:
      print 'Python exception raised in computeF():'
      print e
      return False

  By returning False, you tell NOX that computeF() was unsuccessful.
")
LOCA::Epetra::Interface::Required::computeF;

%import "NOX.Epetra.Interface.i"

%feature("director") LOCA::Epetra::Interface::Required;
%include "LOCA_Epetra_Interface_Required.H"

%feature("director") LOCA::Epetra::Interface::MassMatrix;
%include "LOCA_Epetra_Interface_MassMatrix.H"

%feature("director") LOCA::Epetra::Interface::TimeDependent;
%include "LOCA_Epetra_Interface_TimeDependent.H"

%feature("director") LOCA::Epetra::Interface::TimeDependentMatrixFree;
// The following #define is to change the name of LOCA method
// arguments that conflict with a SWIG director method argument
#define result loca_result
%include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"
#undef result
