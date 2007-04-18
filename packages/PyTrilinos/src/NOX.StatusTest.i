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

    http://software.sandia.gov/trilinos/packages/nox

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
#include "Teuchos_PythonParameter.h"

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
%}

// Trilinos interface file imports
%import "Teuchos.i"

// RefCountPtr typemaps
%teuchos_rcp_typemaps(NOX::StatusTest::Generic)

// General ignore directives
%ignore *::operator=;
%ignore operator<<;
%ignore *::print(ostream& stream, int indent = 0) const;

// Auto-documentation feature
%feature("autodoc", "1");

// SWIG library includes
%include "stl.i"

// NOX::Abstract import
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%import "NOX_Abstract_Group.H"

////////////////////////////////////
// NOX_StatusTest_Generic support //
////////////////////////////////////
%rename(StatusTest_None) NOX::StatusTest::None;
%include "NOX_StatusTest_Generic.H"
namespace NOX {
  namespace StatusTest {
    %extend Generic {
      using namespace std;
      string __str__() {
	stringstream os;
	self->print(os);                  // Put the output in os
	string s = os.str();              // Extract the string from os
	return s.substr(0,s.length()-1);  // Return the string minus trailing \n
      }
    }
  }
}

//////////////////////////////////
// NOX_StatusTest_Combo support //
//////////////////////////////////
%include "NOX_StatusTest_Combo.H"

//////////////////////////////////
// NOX_StatusTest_NormF support //
//////////////////////////////////
%include "NOX_StatusTest_NormF.H"

///////////////////////////////////////
// NOX_StatusTest_NurmUpdate support //
///////////////////////////////////////
%include "NOX_StatusTest_NormUpdate.H"

/////////////////////////////////////
// NOX_StatusTest_NormWRMS support //
/////////////////////////////////////
%include "NOX_StatusTest_NormWRMS.H"

/////////////////////////////////////
// NOX_StatusTest_MaxIters support //
/////////////////////////////////////
%include "NOX_StatusTest_MaxIters.H"

///////////////////////////////////////
// NOX_StatusTest_Stagnation support //
///////////////////////////////////////
%include "NOX_StatusTest_Stagnation.H"

////////////////////////////////////////
// NOX_StatusTest_FiniteValue support //
////////////////////////////////////////
%include "NOX_StatusTest_FiniteValue.H"
