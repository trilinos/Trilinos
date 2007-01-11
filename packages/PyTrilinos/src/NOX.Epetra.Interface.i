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

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1") Interface

// SWIG does not support wrapping nested classes.  We will %import the
// NOX::Utils class (ie, tell swig about it, but not wrap it), which
// has nested classes.  To suppress the swig warning that would
// otherwise result, we use the following:
#pragma SWIG nowarn=312

%{
// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_SrcDistObject.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

// NOX include
#include "NOX_Utils.H"

// NOX::Epetra::Interface includes
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

// Namespace flattening
using namespace NOX;
using namespace NOX::Epetra;
using namespace NOX::Epetra::Interface;
using namespace std;
%}

// Ignore directive
%ignore *::operator=;

// Rename directive
%rename(_print) NOX::Utils::print;

// Feature directives
%feature("director") NOX::Epetra::Interface::Required;
%feature("director") NOX::Epetra::Interface::Jacobian;
%feature("director") NOX::Epetra::Interface::Preconditioner;

// SWIG library includes
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Support for Teuchos::RefCountPtrs
TEUCHOS_RCP_TYPEMAPS(NOX::Epetra::Interface::Required)
TEUCHOS_RCP_TYPEMAPS(NOX::Epetra::Interface::Jacobian)
TEUCHOS_RCP_TYPEMAPS(NOX::Epetra::Interface::Preconditioner)

// NOX imports
%import "NOX_Utils.H"

// NOX::Epetra::Interface includes
%include "NOX_Epetra_Interface_Required.H"
%include "NOX_Epetra_Interface_Jacobian.H"
%include "NOX_Epetra_Interface_Preconditioner.H"
