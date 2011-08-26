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

%define %nox_statustest_docstring
"
PyTrilinos.NOX.StatusTest is the python interface to the namespace
StatusTest in Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.StatusTest is to provide clompletely flexible
specification of stopping tests for NOX algorithms.

NOX.StatusTest provides the following user-level classes:

    * Generic      - Base class for all stopping tests
    * Combo        - Allows combining of stopping tests with AND or OR
    * NormF        - Stopping test based on norm of F
    * NormUpdate   - Stopping test based on norm of update
    * NormWRMS     - Stopping test based on norm of weighted RMS of F
    * MaxIters     - Stopping test based on maximum iterations
    * Stagnation   - Stopping test based on algorithm stagnation
    * FiniteValue  - Stopping test based on detecting NaNs
"
%enddef

%module(package      = "PyTrilinos.NOX",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_statustest_docstring) StatusTest

%{
// System includes
#include <sstream>

// Teuchos includes
#include "PyTrilinos_Teuchos_Util.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include NOX documentation
%include "NOX_dox.i"

// SWIG library includes
%include "stl.i"

// General ignore directives
%ignore *::operator=;
%ignore operator<<;
%ignore *::print(ostream& stream, int indent = 0) const;

// Trilinos interface file imports
%import "Teuchos.i"

// NOX::Abstract import
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%import "NOX_Abstract_Group.H"

// General exception handling
%exception
{
  try
  {
    $action
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unkown C++ exception");
  }
}

// Declare classes to be stored with Teuchos::RCP< >
%teuchos_rcp(NOX::Abstract::Group)
%teuchos_rcp(NOX::Solver::Generic)

////////////////////////////////////
// NOX_StatusTest_Generic support //
////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::Generic)
%rename(StatusTest_None) NOX::StatusTest::None;
%include "NOX_StatusTest_Generic.H"
namespace NOX
{
namespace StatusTest
{
%extend Generic
{
  std::string __str__()
  {
    std::stringstream os;
    self->print(os);                  // Put the output in os
    std::string s = os.str();         // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}
}
}

//////////////////////////////////
// NOX_StatusTest_Combo support //
//////////////////////////////////
%teuchos_rcp(NOX::StatusTest::Combo)
%include "NOX_StatusTest_Combo.H"

//////////////////////////////////
// NOX_StatusTest_NormF support //
//////////////////////////////////
%teuchos_rcp(NOX::StatusTest::NormF)
%include "NOX_StatusTest_NormF.H"

///////////////////////////////////////
// NOX_StatusTest_NurmUpdate support //
///////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::NormUpdate)
%include "NOX_StatusTest_NormUpdate.H"

/////////////////////////////////////
// NOX_StatusTest_NormWRMS support //
/////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::NormWRMS)
%include "NOX_StatusTest_NormWRMS.H"

/////////////////////////////////////
// NOX_StatusTest_MaxIters support //
/////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::MaxIters)
%include "NOX_StatusTest_MaxIters.H"

///////////////////////////////////////
// NOX_StatusTest_Stagnation support //
///////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::Stagnation)
%include "NOX_StatusTest_Stagnation.H"

////////////////////////////////////////
// NOX_StatusTest_FiniteValue support //
////////////////////////////////////////
%teuchos_rcp(NOX::StatusTest::FiniteValue)
%include "NOX_StatusTest_FiniteValue.H"

// Turn off the exception handling
%exception;
