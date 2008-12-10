// -*- c++ -*-

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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %loca_epetra_interface_docstring
"
PyTrilinos.LOCA.Epetra.Interface is the python interface to namespace Epetra::Interface for
the Trilinos package LOCA:

    http://trilinos.sandia.gov/packages/nox

"
%enddef
%module(package      = "PyTrilinos.LOCA.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_interface_docstring) Interface

%{
// NumPy includes
#include "NumPyImporter.h"

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Local Epetra includes
#include "Epetra_NumPyVector.h"

// NOX include
#include "NOX_Utils.H"
#include "NOX_Epetra_Interface_Required.H"

// LOCA::Epetra::Interface includes
#include "LOCA_Epetra_Interface_Required.H"
#include "LOCA_Epetra_Interface_MassMatrix.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"
%}

// General ignore directives
%ignore *::operator=;

// STL support
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCPs typemaps
%teuchos_rcp_typemaps(LOCA::Epetra::Interface::Required)
%teuchos_rcp_typemaps(LOCA::Epetra::Interface::MassMatrix)
%teuchos_rcp_typemaps(LOCA::Epetra::Interface::TimeDependent)
%teuchos_rcp_typemaps(LOCA::Epetra::Interface::TimeDependentMatrixFree)

// Epetra_Vector directorin typemap
%typemap(directorin) Epetra_Vector &
%{
  Epetra_NumPyVector npa$argnum = Epetra_NumPyVector(View,$1_name);
  $input = SWIG_NewPointerObj(&npa$argnum, $descriptor(Epetra_NumPyVector*), 0);
%}

%import "NOX.Epetra.Interface.i"
///////////////////////
// NOX_Utils support //
///////////////////////
// The following #pragma is for nested classes in NOX_Utils.H
#pragma SWIG nowarn=312
%rename(_print) NOX::Utils::print;
%import "NOX_Utils.H"


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
%feature("director") LOCA::Epetra::Interface::Required;
%pythoncode
%{
from NOX.Epetra.Interface import Required
%}
%include "LOCA_Epetra_Interface_Required.H"

%feature("director") LOCA::Epetra::Interface::MassMatrix;
%include "LOCA_Epetra_Interface_MassMatrix.H"

%feature("director") LOCA::Epetra::Interface::TimeDependent;
%include "LOCA_Epetra_Interface_TimeDependent.H"

%feature("director") LOCA::Epetra::Interface::TimeDependentMatrixFree;
// The following #define is to change the name of LOCA method
// arguments that conflict with a SWIG director method argument
#define result nox_result
%include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"

