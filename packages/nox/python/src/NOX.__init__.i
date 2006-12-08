// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package="PyTrilinos.NOX") __init__

%{
// System includes
#include <sstream>

// Teuchos include
#include "Teuchos_PythonParameter.hpp"

// NOX includes
#include "NOX_Version.H"
#include "NOX_Utils.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_Manager.H"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"
%}

// Ignore directives
%ignore operator<<;
%ignore *::operator=;
%ignore NOX::Utils::fill;
%ignore NOX::Utils::sciformat;

// Rename directive
%rename(_print) NOX::Utils::print;

// Auto-documentation feature
%feature("autodoc", "1");

// SWIG library includes
%include "stl.i"

// Trilinos interface file imports.  Note: Teuchos.i turns off
// warnings for nested classes, so we do not have to do it again.
%import "Teuchos.i"

// NOX top-level interface includes
using namespace std;
%include "NOX_Version.H"
%include "NOX_Utils.H"

// NOX namespace imports
%import "NOX.Abstract.i"
%import "NOX.Solver.i"
%import "NOX.StatusTest.i"

// Python code for the NOX module
%pythoncode %{
__version__ = version().split()[2]

try:
    import Epetra
except ImportError:
    pass
%}
