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

%define %loca_docstring
"
PyTrilinos.LOCA is the python interface to the Trilinos continuation
algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA is to provide a library of continuation
algorithms.  This module is not currently supported, but the plan is
to reactivate it soon.
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring = %loca_docstring) __init__

%{
// System includes
#include <sstream>

// Teuchos include
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"

// LOCA includes
#include "LOCA.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"

//#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
//#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
//#include "LOCA_MultiContinuation_AbstractGroup.H"
//#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

#include "LOCA_TimeDependent_AbstractGroup.H"
#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"

#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

// Local includes
#include "NumPyImporter.h"

// Namespace flattening
using Teuchos::RCP;

%}

// Ignore/renames
%ignore *::operator=;
%ignore *::operator[];
%ignore operator<<(ostream&, const LOCA::ParameterVector&);
%rename(Print) LOCA::ParameterVector::print(ostream& stream) const;

// SWIG library includes
%include "stl.i"

// Trilinos interface file imports.
// The Teuchos.py file that this %import will try to import in python
// will be one directory up from this python module.  So we need to
// add the parent directory to the search path.
%pythoncode
{
import os.path, sys
currentDir,dummy = os.path.split(__file__)
sys.path.append(os.path.normpath(os.path.join(currentDir,"..")))
}
%import "Teuchos.i"
// Note: Teuchos.i turns off warnings for nested classes, so we do not
// have to do it again.

%import "LOCA.MultiContinuation.i"
// %import "LOCA.Continuation.i"

//%import "NOX.StatusTest.i"
//%include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
//%include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"

// NOX interface file imports.
%pythoncode
{
import os.path, sys
currentDir,dummy = os.path.split(__file__)
sys.path.append(os.path.normpath(os.path.join(currentDir,"..","NOX")))
}
//%import "NOX.__init__.i"

//%include "LOCA_MultiContinuation_AbstractGroup.H"
//%include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

%rename(TimeDependent_AbstractGroup) LOCA::TimeDependent::AbstractGroup;
%include "LOCA_TimeDependent_AbstractGroup.H"
%rename(Homotopy_AbstractGroup) LOCA::Homotopy::AbstractGroup;
%include "LOCA_Homotopy_AbstractGroup.H"
%rename(TurningPoint_MooreSpence_AbstractGroup) LOCA::TurningPoint::MooreSpence::AbstractGroup;
%include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%rename(TurningPoint_MinimallyAugmented_AbstractGroup) LOCA::TurningPoint::MinimallyAugmented::AbstractGroup;
%include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
%rename(TurningPoint_MooreSpence_FiniteDifferenceGroup) LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup;
%include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%rename(TurningPoint_MinimallyAugmented_FiniteDifferenceGroup) LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup;
%include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
%rename(Pitchfork_MooreSpence_AbstractGroup) LOCA::Pitchfork::MooreSpence::AbstractGroup;
%include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%rename(Pitchfork_MinimallyAugmented_AbstractGroup) LOCA::Pitchfork::MinimallyAugmented::AbstractGroup;
%include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
%rename(Hopf_MooreSpence_AbstractGroup) LOCA::Hopf::MooreSpence::AbstractGroup;
%include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%rename(Hopf_MinimallyAugmented_AbstractGroup) LOCA::Hopf::MinimallyAugmented::AbstractGroup;
%include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
%rename(Hopf_MooreSpence_FiniteDifferenceGroup) LOCA::Hopf::MooreSpence::FiniteDifferenceGroup;
%include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
%rename(Hopf_MinimallyAugmented_FiniteDifferenceGroup) LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup;
%include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"


//%rename(Abstract_Group) LOCA::Abstract::Group;
//%include "LOCA_Abstract_Group.H"
//%rename(Abstract_TransposeSolveGroup) LOCA::Abstract::TransposeSolveGroup;
//%include "LOCA_Abstract_TransposeSolveGroup.H"
//%include "LOCA_Abstract_Iterator.H"

// LOCA interface includes
%include "LOCA.H"
%include "LOCA_GlobalData.H"
%teuchos_rcp_typemaps(LOCA::GlobalData)
%teuchos_rcp_typemaps(LOCA::DerivUtils)

%import "LOCA_Abstract_Iterator.H"
%import "NOX.StatusTest.i"

%include "LOCA_Stepper.H"
%include "LOCA_Parameter_Vector.H"


%pythoncode
%{
import Epetra
%}
