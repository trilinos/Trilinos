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

%define %loca_abstract_docstring
"
PyTrilinos.LOCA.Abstract is the python interface to namespace Abstract
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Abstract is to provide abstract continuation
problem base classes.  The python version of LOCA.Abstract supports
the following classes:

    * Group                - Compatiblity class for AbstractGroup hierarchy
    * TransposeSolveGroup  - Abstract group interface class for solving the
                             transpose of the Jacobian
    * Iterator             - Abstract interface for implementing iteration
    * Factory              - Abstract interface for providing a user-defined
                             factory
"
%enddef

%module(package      = "PyTrilinos.LOCA",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_abstract_docstring) Abstract

%{
// PyTrilinos includes
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_Teuchos_Util.hpp"

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Director exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}

// General ignore directives
%ignore *::operator=;
%ignore operator=;
%ignore *::operator[];

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Homotopy::AbstractGroup)
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Pitchfork::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Pitchfork::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::PhaseTransition::AbstractGroup)
%teuchos_rcp(LOCA::Abstract::Group)
%teuchos_rcp(LOCA::Abstract::TransposeSolveGroup)
%teuchos_rcp(LOCA::Abstract::Iterator)
%teuchos_rcp(LOCA::Abstract::Factory)

// Import SWIG interface files to provide information about base
// classes
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="Homotopy") "LOCA_Homotopy_AbstractGroup.H"
%import(module="TimeDependent") "LOCA_TimeDependent_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="Pitchfork.MooreSpence") "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%import(module="Pitchfork.MinimallyAugmented") "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="PhaseTransition") "LOCA_PhaseTransition_AbstractGroup.H"

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
  SWIG_CATCH_STDEXCEPT
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch (Swig::DirectorMethodException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// LOCA::Abstract Group class
// The following #define is to change the name of LOCA method
// arguments that conflict with a SWIG director method argument
#define result loca_result
%include "LOCA_Abstract_Group.H"

// LOCA::Abstract TransposeSolveGroup class
%include "LOCA_Abstract_TransposeSolveGroup.H"

// LOCA::Abstract Iterator class
%include "LOCA_Abstract_Iterator.H"

// LOCA::Abstract Factory class
%include "LOCA_Abstract_Factory.H"

// We need to undefine the result macro
#undef result
