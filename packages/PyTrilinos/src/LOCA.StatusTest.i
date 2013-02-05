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

%module(package="PyTrilinos.LOCA") StatusTest

%{
// PyTrilinos include
#include "PyTrilinos_Teuchos_Util.h"

// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// // Configuration and optional includes
// %include "PyTrilinos_config.h"
// #ifdef HAVE_NOX_EPETRA
// %{
// #include "NOX_Epetra_Group.H"
// #include "NOX_Epetra_Vector.H"
// #include "Epetra_NumPyVector.h"
// %}
// #endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"
%import "NOX.StatusTest.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::StatusTest::Abstract)
%teuchos_rcp(LOCA::StatusTest::Wrapper)

// LOCA::StatusTest Abstract class
%rename (None_) None;
%include "LOCA_StatusTest_Abstract.H"

// LOCA::StatusTest Wrapper class
%include "LOCA_StatusTest_Wrapper.H"
