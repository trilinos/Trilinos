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
#include "PyTrilinos_Teuchos_Util.h"

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
#include "LOCA_Abstract_Iterator.H"
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"

//#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
//#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
//#include "LOCA_MultiContinuation_AbstractGroup.H"
//#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"

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
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

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

// Trilinos interface import
%import "Teuchos.i"
// Note: Teuchos.i turns off warnings for nested classes, so we do not
// have to do it again.

// Exception handling
%include "exception.i"

// Director exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PyTrilinos::PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  catch(int errCode)
  {
    PyErr_Format(PyExc_EpetraError, "Error code = %d\nSee stderr for details", errCode);
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

%teuchos_rcp(LOCA::GlobalData)
%teuchos_rcp(LOCA::DerivUtils)

// NOX interface file imports.
//%import "NOX.__init__.i"
%import "NOX.Abstract.i"
%import "NOX.StatusTest.i"

%import "LOCA.Abstract.i"
%import "LOCA.Extended.i"
%import "LOCA.BorderedSystem.i"
%import "LOCA.Continuation.i"
%import "LOCA.MultiContinuation.i"
%import "LOCA.Hopf.i"
%import "LOCA.TimeDependent.i"
%import "LOCA.Pitchfork.i"
%import "LOCA.Homotopy.i"
%import "LOCA.TurningPoint.i"
//%import "LOCA_Abstract_Iterator.H"

// LOCA interface includes
%include "LOCA.H"
%include "LOCA_GlobalData.H"


%include "LOCA_Stepper.H"
%include "LOCA_Parameter_Vector.H"


%pythoncode
%{
#import Epetra
%}
