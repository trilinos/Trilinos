// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
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
// System include files
#include <sstream>

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// NOX include files
#include "PyTrilinos_NOX_Abstract_Headers.hpp"
#include "PyTrilinos_NOX_StatusTest_Headers.hpp"
#include "PyTrilinos_NOX_Solver_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include NOX documentation
%include "NOX_dox.i"

// SWIG library include files
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
