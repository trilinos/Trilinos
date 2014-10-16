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

%define %loca_anasazioperator_docstring
"
PyTrilinos.LOCA.AnasaziOperator is the python interface to namespace
AnasaziOperator of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.AnasaziOperator is to provide ***.  The python
version of LOCA.AnasaziOperator supports the following classes:

    * AbstractStrategy  - Abstract interface class for Anasazi operator
                          strategies
    * Factory           - Factory for creating Anasazi operator strategy
                          objects
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_anasazioperator_docstring) AnasaziOperator

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
#include "LOCA_AnasaziOperator_AbstractStrategy.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Namespace flattening
using Teuchos::RCP;
%}

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::AnasaziOperator::AbstractStrategy)
%teuchos_rcp(LOCA::AnasaziOperator::Factory)

// LOCA::AnasaziOperator AbstractStrategy class
%feature("director") LOCA::AnasaziOperator::AbstractStrategy;
%extend LOCA::AnasaziOperator::AbstractStrategy
{
  // The C++ version of this method returns a reference to a
  // std::string.  This is not recommended for director methods, so I
  // re-implement it here to return a copy of the label.
  const std::string label() const
  {
    return self->label();
  }
}
%ignore LOCA::AnasaziOperator::AbstractStrategy::label;
%include "LOCA_AnasaziOperator_AbstractStrategy.H"

// LOCA::AnasaziOperator Factory class
%include "LOCA_AnasaziOperator_Factory.H"
