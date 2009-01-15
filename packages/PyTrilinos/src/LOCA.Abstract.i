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

%define %loca_abstract_docstring
"
PyTrilinos.LOCA.Abstract is the python interface to namespace Abstract
of the Trilinos package LOCA:

    http://trilinos.sandia.gov/packages/nox

"
%enddef

%module(package      = "PyTrilinos.LOCA",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_abstract_docstring) Abstract

%{
// Teuchos includes
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"

// LOCA includes
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_Iterator.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_Stepper.H"

// Local includes
#include "NumPyImporter.h"
%}

// Standard exception handling
%include "exception.i"

// Include NOX documentation
// %include "LOCA_dox.i"   // TODO: this file will need to be generated

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

%pythoncode
{
import os.path, sys
currentDir,dummy = os.path.split(__file__)
sys.path.append(os.path.normpath(os.path.join(currentDir,"..")))
}

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCPs typemaps
//%teuchos_rcp_typemaps(LOCA::GlobalData)
//%teuchos_rcp_typemaps(LOCA::DerivUtils)

%include "LOCA_Abstract_Iterator.H"

%import "NOX.Abstract.i"
%import "LOCA.Homotopy.i"
%import "LOCA.__init__.i"

%include "LOCA_Abstract_Group.H"
%include "LOCA_Abstract_TransposeSolveGroup.H"
