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

%define %nox_solver_docstring
"
PyTrilinos.NOX.Solver is the python interface to the Solver namespace
of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Solver is to provide solver manager classes for
NOX.  NOX.Solver provides the following user-level classes:

    * Generic                  - Base class for solver managers
    * LineSearchBased          - Line-search-based solver manager
    * TrustRegionBased         - Trust-region-based solver manager
    * InexactTrustRegionBased  - Inexact-trust-region-based solver
                                 manager
    * TensorBased              - Tensor-based solver manager

in addition to the following factory function:

    * buildSolver              - Recommended method for creating solver
                                 managers (note that without loss of
                                 functionality, the Factory class is not
                                 currently provided).
"
%enddef

%module(package      = "PyTrilinos.NOX",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_solver_docstring) Solver

%{
// Teuchos includes
#include "PyTrilinos_Teuchos_Util.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"
#include "NOX_Solver_Factory.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

// Namespace flattening
using Teuchos::RCP;
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
%}
#endif

// Standard exception handling
%include "exception.i"

// Include NOX documentation
%include "NOX_dox.i"

// General ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream& stream, int indent = 0) const;

// Rename directives
%rename(StatusTest_Generic) NOX::StatusTest::Generic;
%rename(StatusTest_None   ) NOX::StatusTest::None;

// Trilinos imports
%import "Teuchos.i"
%import "NOX.Abstract.i"
%import "NOX.StatusTest.i"

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType & e)
  {
    SWIG_exception(SWIG_TypeError, e.what());
  }
  catch(Teuchos::Exceptions::InvalidParameter & e)
  {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  catch(const char* e)
  {
    PyErr_SetString(PyExc_RuntimeError, e);
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Downcast NOX::Abstract::Group return arguments to NOX::Epetra::Group,
// if possible
#ifdef HAVE_NOX_EPETRA
%typemap(out) const NOX::Abstract::Group & (NOX::Epetra::Group* negResult = NULL)
{
  negResult = dynamic_cast< NOX::Epetra::Group* >(const_cast< NOX::Abstract::Group* >($1));
  if (negResult)
  {
    Teuchos::RCP< NOX::Epetra::Group > *smartresult = new
      Teuchos::RCP< NOX::Epetra::Group >(negResult, bool($owner));
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Epetra::Group > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    // If we cannot downcast, then return the NOX::Abstract::Group
    Teuchos::RCP< NOX::Abstract::Group > *smartresult =
      new Teuchos::RCP< NOX::Abstract::Group >($1, bool($owner));
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Abstract::Group > * ),
				   SWIG_POINTER_OWN));
  }
}
#endif

//////////////////////////////////
// NOX::Solver::Generic support //
//////////////////////////////////
%ignore *::getSolutionGroup;
%ignore *::getPreviousSolutionGroup;
%ignore *::getList;
%rename(getSolutionGroup        ) *::getSolutionGroupPtr;
%rename(getPreviousSolutionGroup) *::getPreviousSolutionGroupPtr;
%rename(getList                 ) *::getListPtr;
%teuchos_rcp(NOX::Solver::Generic)
%include "NOX_Solver_Generic.H"

//////////////////////////////////////////
// NOX::Solver::LineSearchBased support //
//////////////////////////////////////////
%teuchos_rcp(NOX::Solver::LineSearchBased)
%include "NOX_Solver_LineSearchBased.H"

///////////////////////////////////////////
// NOX::Solver::TrustRegionBased support //
///////////////////////////////////////////
%teuchos_rcp(NOX::Solver::TrustRegionBased)
%include "NOX_Solver_TrustRegionBased.H"

//////////////////////////////////////////////////
// NOX::Solver::InexactTrustRegionBased support //
//////////////////////////////////////////////////
%teuchos_rcp(NOX::Solver::InexactTrustRegionBased)
%include "NOX_Solver_InexactTrustRegionBased.H"

//////////////////////////////////////
// NOX::Solver::TensorBased support //
//////////////////////////////////////
%teuchos_rcp(NOX::Solver::TensorBased)
%include "NOX_Solver_TensorBased.H"

//////////////////////////////////////
// NOX::Solver::buildSolver support //
//////////////////////////////////////
%rename (buildSolver) myBuildSolver;
// NOX::Solver::buildSolver in NOX_Solver_Factory.H returns a
// Teuchos::RCP<NOX::Solver::Generic>.  As far as I can tell, SWIG
// cannot properly downcast the NOX::Solver::Generic object wrapped
// within the Teuchos::RCP<> in order to, say, call its solve()
// method.  Therefore, I write my own wrapper around buildSolver()
// that does this downcasting explicitly and returns a python wrapper
// around the appropriate derived class.
%inline
{
  PyObject *
    myBuildSolver(const Teuchos::RCP< NOX::Abstract::Group     > & grp,
		  const Teuchos::RCP< NOX::StatusTest::Generic > & tests,
		  const Teuchos::RCP< Teuchos::ParameterList   > & params)
  {
    // SWIG type queries
    static swig_type_info * swig_NSLSB_ptr  =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::LineSearchBased > *");
    static swig_type_info * swig_NSTRB_ptr  =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::TrustRegionBased > *");
    static swig_type_info * swig_NSITRB_ptr =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > *");
    static swig_type_info * swig_NSTB_ptr   =
      SWIG_TypeQuery("Teuchos::RCP< NOX::Solver::TensorBased > *");
    // Build a NOX::Solver::Generic object via the buildSolver factory
    Teuchos::RCP< NOX::Solver::Generic > rcp_solver = NOX::Solver::buildSolver(grp, tests, params);
    // Try to downcast to a derived class
    {
      Teuchos::RCP< NOX::Solver::LineSearchBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::LineSearchBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::LineSearchBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::LineSearchBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSLSB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::TrustRegionBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::TrustRegionBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::TrustRegionBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::TrustRegionBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSTRB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::InexactTrustRegionBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::InexactTrustRegionBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::InexactTrustRegionBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSITRB_ptr, 1);
      }
    }
    {
      Teuchos::RCP< NOX::Solver::TensorBased > result =
	Teuchos::rcp_dynamic_cast< NOX::Solver::TensorBased >(rcp_solver);
      if (!result.is_null())
      {
	Teuchos::RCP< NOX::Solver::TensorBased > *smartresult =
	  new Teuchos::RCP< NOX::Solver::TensorBased >(result);
	return SWIG_NewPointerObj((void*)smartresult, swig_NSTB_ptr, 1);
      }
    }
    PyErr_SetString(PyExc_RuntimeError, "NOX::Solver::buildSolver returned unrecognized "
		    "derivative of NOX::Solver::Generic");
    return NULL;
  }
}

// Turn off the exception handling
%exception;
