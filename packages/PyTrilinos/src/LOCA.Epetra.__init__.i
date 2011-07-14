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

%define %loca_epetra_docstring
"
PyTrilinos.LOCA.Epetra is the python interface to namespace Epetra for
the Trilinos package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra is to provide a concrete interface beteen
LOCA and Epetra.

"
%enddef
%module(package      = "PyTrilinos.LOCA.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_docstring) __init__

%{
// System includes
#include <vector>

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FEVbrMatrix.h"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"

//#include "NOX_Abstract_Group.H"
#include "NOX_Epetra_Group.H"
//#include "NOX_Epetra_Interface_Preconditioner.H"
//#include "NOX_Epetra_FiniteDifference.H"
//#include "NOX_Epetra_FiniteDifferenceColoring.H"
//#include "NOX_Epetra_LinearSystem_AztecOO.H"
//#include "NOX_Epetra_MatrixFree.H"
//#include "LOCA_MultiContinuation_AbstractGroup.H"
//#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_TimeDependent_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"
#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"
#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_BorderedSystem_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"

#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
#include "LOCA_Epetra.H"
#include "LOCA_Epetra_Group.H"

// Namespace flattening
using Teuchos::RCP;
using Teuchos::rcp;
%}

%ignore *::operator=;

// SWIG library includes
%include "stl.i"

// Trilinos interface support
%import "Teuchos.i"

%import "NOX.Abstract.i"
%import "NOX.Epetra.__init__.i"

%import "LOCA.__init__.i"
%import "LOCA.MultiContinuation.i"

//%import "LOCA_TimeDependent_AbstractGroup.H"
//%import "LOCA_Homotopy_AbstractGroup.H"
//%import "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
//%import "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
//%import "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
//%import "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
//%import "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
//%import "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
//%import "LOCA_Hopf_MooreSpence_AbstractGroup.H"
//%import "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
//%import "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
//%import "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"

//%import "LOCA.Abstract.i"
%import "LOCA.Hopf.i"
%import "LOCA.Pitchfork.i"
%import "LOCA.Homotopy.i"
%import "LOCA.TurningPoint.i"
// %pythoncode
// %{
// import PyTrilinos.LOCA.Homotopy
// import PyTrilinos.LOCA.TurningPoint
// import PyTrilinos.LOCA.Pitchfork
// import PyTrilinos.LOCA.Hopf
// %}

%import "NOX.Epetra.Interface.i"
%import "LOCA.Epetra.Interface.i"

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

%rename(Abstract_Group) LOCA::Abstract::Group;
%rename(Abstact_TransposeSolveGroup) LOCA::Abstract::TransposeSolveGroup;
%teuchos_rcp(Abstract_Group)
%teuchos_rcp(Abstract_TransposeSolveGroup)
%include "LOCA_Abstract_Group.H"
%include "LOCA_Abstract_TransposeSolveGroup.H"

#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
%include "LOCA_Epetra.H"

// %pythoncode
// %{
// from NOX.Epetra import Group
// %}

//////////////////////////////
// LOCA.Epetra.Group support //
//////////////////////////////

%rename(NOX_Epetra_Group) NOX::Epetra::Group;
%include "NOX_Epetra_Group.H"

// temporarily ignore conflict-causing constructor.  TODO: fix this issue
%ignore LOCA::Epetra::Group::Group(Teuchos::RCP< LOCA::GlobalData > const &,Teuchos::ParameterList &,Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree > const &,NOX::Epetra::Vector &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,LOCA::ParameterVector const &);

%include "LOCA_Epetra_Group.H"

%teuchos_rcp(LOCA::Epetra::Group)
