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

%define %nox_epetra_interface_docstring
"
PyTrilinos.NOX.Epetra.Interface is the python interface to the
Epetra::Interface namespace of the Trilinos package NOX:

    http://software.sandia.gov/trilinos/packages/nox

The purpose of NOX.Epetra.Interface is to provide base classes the
user should derive from in order to define the nonlinear function to
be solved, and if needed, its Jacobian and the desired preconditioner.

NOX.Epetra.Interface provides the following user-level classes:

    * Required        - Required class for computing the nonlinear function
    * Jacobian        - Class for computing the Jacobian (if needed)
    * Preconditioner  - Class for computing the preconditioner (if needed)
"
%enddef

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_epetra_interface_docstring) Interface

%{
// NumPy includes
#include "NumPyImporter.h"

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_SrcDistObject.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

// Local Epetra includes
#include "Epetra_NumPyVector.h"

// NOX include
#include "NOX_Utils.H"

// NOX::Epetra::Interface includes
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
%}

// General ignore directives
%ignore *::operator=;

// STL support
using namespace std;
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RefCountPtrs typemaps
%teuchos_rcp_typemaps(NOX::Epetra::Interface::Required)
%teuchos_rcp_typemaps(NOX::Epetra::Interface::Jacobian)
%teuchos_rcp_typemaps(NOX::Epetra::Interface::Preconditioner)

// Epetra_Vector directorin typemap
%typemap(directorin) Epetra_Vector & %{
  Epetra_NumPyVector npa$argnum = Epetra_NumPyVector(View,$1_name);
  $input = SWIG_NewPointerObj(&npa$argnum, $descriptor(Epetra_NumPyVector*), 0);
%}

///////////////////////
// NOX_Utils support //
///////////////////////
// The following #pragma is for nested classes in NOX_Utils.H
#pragma SWIG nowarn=312
%rename(_print) NOX::Utils::print;
%import "NOX_Utils.H"

///////////////////////////////////////////
// NOX_Epetra_Interface_Required support //
///////////////////////////////////////////
%feature("director") NOX::Epetra::Interface::Required;
%include "NOX_Epetra_Interface_Required.H"

///////////////////////////////////////////
// NOX_Epetra_Interface_Jacobian support //
///////////////////////////////////////////
%feature("director") NOX::Epetra::Interface::Jacobian;
%include "NOX_Epetra_Interface_Jacobian.H"

/////////////////////////////////////////////////
// NOX_Epetra_Interface_Preconditioner support //
/////////////////////////////////////////////////
%feature("director") NOX::Epetra::Interface::Preconditioner;
%include "NOX_Epetra_Interface_Preconditioner.H"
