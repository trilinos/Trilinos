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

%define %loca_hopf_docstring
"
PyTrilinos.LOCA.Hopf is the python interface to namespace Hopf of the
Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Hopf is to provide groups and vectors for Hopf
bifurcations.  The python version of LOCA.Hopf supports the following
sub-modules:

    * MooreSpence         - Groups and vectors for locating Hopf bifurcations
                            using the Moore-Spence formulation
    * MinimallyAugmented  - Groups and vectors for locating Hopf bifurcations
                            using the minimally augmented Hopf formulation

and classes:

    * ComplexMultiVector  - Multi-vector class to hold two multi-vectors to
                            represent a complex multi-vector
    * ComplexVector       - Vector class to hold two vectors to represent a
                            complex vector
"
%enddef

%module(package   = "PyTrilinos.LOCA.Hopf",
        docstring = %loca_hopf_docstring) __init__

%{
// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_ComplexMultiVector.H"
#include "LOCA_Hopf_ComplexVector.H"

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
%ignore *::operator[];
%ignore operator[];

// Trilinos module imports
%import "Teuchos.i"

// Base class imports
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}
%import "NOX.Abstract.i"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"

// Import the sub-modules
%pythoncode
%{
__all__ = ['MooreSpence',
           'MinimallyAugmented'
           ]
import MooreSpence
import MinimallyAugmented
%}

// LOCA::Hopf ComplexMultiVector class
%include "LOCA_Hopf_ComplexMultiVector.H"

// LOCA::Hopf ComplexVector class
%include "LOCA_Hopf_ComplexVector.H"
