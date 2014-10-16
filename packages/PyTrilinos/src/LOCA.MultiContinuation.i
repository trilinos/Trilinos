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

%define %loca_multicontinuation_docstring
"
PyTrilinos.LOCA.MultiContinuation is the python interface to namespace
MultiContinuation of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.MultiContinuation is to provide groups and vectors
for multi-parameter continuation.  The python version of
LOCA.MultiContinuation supports the following classes:

    * AbstractGroup            - LOCA abstract interface for continuation,
                                 derived from the NOX.Abstract.Group.  This
                                 abstract class provides the interface
                                 necessary to perform continuation, i.e.,
                                 compute families of solutions to F(x,p) = 0
    * FiniteDifferenceGroup    - Concrete class that provides a concrete
                                 implementation of the computeDfDp() method of
                                 the LOCA.Continuation.AbstractGroup using
                                 first-order finite differencing
    * ConstraintInterface      - Abstract interface for the constraint portion
                                 of a constrained nonlinear system
    * ConstraintInterfaceMVDX  - Abstract interface for the constraint portion
                                 of a constrained nonlinear system for
                                 constraints that support computing a solution
                                 component derivative as a multi-vector
    * ExtendedMultiVector      - MultiVector class to hold solution vectors,
                                 Newton vectors, etc. for continuation equations
    * ExtendedVector           - Vector class to hold solution vectors, Newton
                                 vectors, etc. for continuation equations
    * Factory                  - Factory for creating continuation strategy
                                 objects
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_multicontinuation_docstring) MultiContinuation

%{
// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_NumPyVector.hpp"
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
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::MultiContinuation::Factory)

// Import base class declarations
%import "NOX.Abstract.i"
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"

// LOCA::MultiContinuation AbstractGroup class
%include "LOCA_MultiContinuation_AbstractGroup.H"

// LOCA::MultiContinuation FiniteDifferenceGroup class
%include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

// LOCA::MultiContinuation ConstraintInterface class
%feature("director") LOCA::MultiContinuation::ConstraintInterface;
%include "LOCA_MultiContinuation_ConstraintInterface.H"

// LOCA::MultiContinuation ConstraintInterfaceMVDX class
%warnfilter(473) LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%feature("director") LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

// LOCA::MultiContinuation ExtendedMultiVector class
%include "LOCA_MultiContinuation_ExtendedMultiVector.H"

// LOCA::MultiContinuation ExtendedVector class
%include "LOCA_MultiContinuation_ExtendedVector.H"

// LOCA::MultiContinuation Factory class
%include "LOCA_MultiContinuation_Factory.H"
