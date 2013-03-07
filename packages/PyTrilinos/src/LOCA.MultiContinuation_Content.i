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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
// Teuchos includes
#include "PyTrilinos_Teuchos_Util.h"

// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_NumPyVector.h"
%}
#endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

%import "Teuchos.i"

// Teucho::RCP support
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::MultiContinuation::Factory)

// Import base class declarations
%import "NOX.Abstract.i"
%import "LOCA.Extended.i"

// LOCA::MultiContinuation AbstractGroup class
%feature("director") LOCA::MultiContinuation::AbstractGroup;
%include "LOCA_MultiContinuation_AbstractGroup.H"

// LOCA::MultiContinuation FiniteDifferenceGroup class
%feature("director") LOCA::MultiContinuation::FiniteDifferenceGroup;
%include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

// LOCA::MultiContinuation ConstraintInterface class
%feature("director") LOCA::MultiContinuation::ConstraintInterface;
%include "LOCA_MultiContinuation_ConstraintInterface.H"

// LOCA::MultiContinuation ConstraintInterfaceMVDX class
%feature("director") LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

// LOCA::MultiContinuation ExtendedMultiVector class
%feature("director") LOCA::MultiContinuation::ExtendedMultiVector;
%include "LOCA_MultiContinuation_ExtendedMultiVector.H"

// LOCA::MultiContinuation ExtendedVector class
%feature("director") LOCA::MultiContinuation::ExtendedVector;
%include "LOCA_MultiContinuation_ExtendedVector.H"

// LOCA::MultiContinuation Factory class
%include "LOCA_MultiContinuation_Factory.H"
