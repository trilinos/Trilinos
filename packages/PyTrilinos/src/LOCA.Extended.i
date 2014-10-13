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

%define %loca_extended_docstring
"
PyTrilinos.LOCA.Extended is the python interface to namespace Extended
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Extended is to provide classes that extend
NOX.Abstract classes to handle an arbitrary number of multi-vectors
and scalars.  The python version of LOCA.Extended supports the
following classes:

    * MultiVector         - Implemenatation of the NOX.Abstract.MultiVector
                            class for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * Vector              - Implemenatation of the NOX.Abstract.Vector class
                            for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * MultiAbstractGroup  - LOCA abstract interface for extended groups,
                            derived from the NOX.Abstract.Group, i.e., an
                            abstract interface for 'super' groups that have an
                            underlying group component
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_extended_docstring) Extended

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

// Trilinos module imports
%import "Teuchos.i"
%import "NOX.Abstract.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)

//////////////////////////////////////
// LOCA::Extended MultiVector class //
//////////////////////////////////////
//%feature("director") LOCA::Extended::MultiVector;
%include "LOCA_Extended_MultiVector.H"

/////////////////////////////////
// LOCA::Extended Vector class //
/////////////////////////////////
//%feature("director") LOCA::Extended::Vector;
%ignore LOCA::Extended::Vector::getScalar(int);
%include "LOCA_Extended_Vector.H"

/////////////////////////////////////////////
// LOCA::Extended MultiAbstractGroup class //
/////////////////////////////////////////////
//%feature("director") LOCA::Extended::MultiAbstractGroup;
%include "LOCA_Extended_MultiAbstractGroup.H"
