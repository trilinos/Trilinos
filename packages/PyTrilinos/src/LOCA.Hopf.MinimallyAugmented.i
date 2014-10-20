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

%define %loca_hopf_minimallyaugmented_docstring
"
PyTrilinos.LOCA.Hopf.MinimallyAugmented is the python interface to
namespace Hopf::MinimallyAugmented of the Trilinos continuation
algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Hopf.MinimallyAugmented is to provide groups and
vectors for locating Hopf bifurcations using the minimally augmented
Hopf formulation.  The python version of LOCA.Hopf.MinimallyAugmented
supports the following classes:

    * Constraint             - Implementation of LOCA.MultiContinuation.-
                               ConstraintInterfaceMVDX for computing Hopf
                               bifurcations for the minimally augmented Hopf
                               formulation
    * AbstractGroup          - Interface to underlying groups for Hopf
                               calculations using the minimally augmented
                               formulation
    * FiniteDifferenceGroup  - Concrete class that provides concrete
                               implementations of the derivative computation
                               methods of the LOCA.Hopf.MinimallyAugmented.-
                               AbstractGroup using first-order finite
                               differencing
    * ExtendedGroup          - Group representing the minimally augemented Hopf
                               equations
"
%enddef

%module(package   = "PyTrilinos.LOCA.Hopf",
        directors = "1",
        docstring = %loca_hopf_minimallyaugmented_docstring) MinimallyAugmented

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
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"

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

// Trilinos module imports
%import "Teuchos.i"

// The following #define is to change the name of LOCA method
// arguments that conflict with a SWIG director method argument
#define result loca_result

// Teuchos::RCP handling
%teuchos_rcp(LOCA::BorderedSystem::AbstractGroup)
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::Constraint)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::ExtendedGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)

// Base class support
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}
%import "NOX.Abstract.i"
%import(module="BorderedSystem") "LOCA_BorderedSystem_AbstractGroup.H"
%warnfilter(473) LOCA::MultiContinuation::AbstractGroup;
%warnfilter(473) LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterface.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
%import(module="TimeDependent") "LOCA_TimeDependent_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%import(module="MooreSpence") "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%import(module="MooreSpence") "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::Hopf::MinimallyAugmented Constraint class
%warnfilter(473) LOCA::Hopf::MinimallyAugmented::Constraint;
%feature("director") LOCA::Hopf::MinimallyAugmented::Constraint;
%include "LOCA_Hopf_MinimallyAugmented_Constraint.H"

// LOCA::Hopf::MinimallyAugmented AbstractGroup class
%feature("director") LOCA::Hopf::MinimallyAugmented::AbstractGroup;
%include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"

// LOCA::Hopf::MinimallyAugmented FiniteDifferenceGroup class
%feature("director") LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup;
%include "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"

// LOCA::Hopf::MinimallyAugmented ExtendedGroup class
%feature("director") LOCA::MinimallyAugmented::ExtendedGroup;
%include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"

// We need to clear this
#undef loca_result
