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
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

//////////////////////////////////
// NOX::Solver::Generic support //
//////////////////////////////////
%include "NOX_Solver_Generic.H"
%teuchos_rcp_typemaps(NOX::Solver::Generic)

//////////////////////////////////////////
// NOX::Solver::LineSearchBased support //
//////////////////////////////////////////
%include "NOX_Solver_LineSearchBased.H"
%teuchos_rcp_typemaps(NOX::Solver::LineSearchBased)

///////////////////////////////////////////
// NOX::Solver::TrustRegionBased support //
///////////////////////////////////////////
%include "NOX_Solver_TrustRegionBased.H"
%teuchos_rcp_typemaps(NOX::Solver::TrustRegionBased)

//////////////////////////////////////////////////
// NOX::Solver::InexactTrustRegionBased support //
//////////////////////////////////////////////////
%include "NOX_Solver_InexactTrustRegionBased.H"
%teuchos_rcp_typemaps(NOX::Solver::InexactTrustRegionBased)

//////////////////////////////////////
// NOX::Solver::TensorBased support //
//////////////////////////////////////
%include "NOX_Solver_TensorBased.H"
%teuchos_rcp_typemaps(NOX::Solver::TensorBased)

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
    myBuildSolver(const Teuchos::RCP<NOX::Abstract::Group>     & grp,
		  const Teuchos::RCP<NOX::StatusTest::Generic> & tests,
		  const Teuchos::RCP<Teuchos::ParameterList>   & params)
  {
    // SWIG type queries
    static swig_type_info * swig_NSLSB_ptr  =
      SWIG_TypeQuery("NOX::Solver::LineSearchBased*");
    static swig_type_info * swig_NSTRB_ptr  =
      SWIG_TypeQuery("NOX::Solver::TrustRegionBased*");
    static swig_type_info * swig_NSITRB_ptr =
      SWIG_TypeQuery("NOX::Solver::InexactTrustRegionBased*");
    static swig_type_info * swig_NSTB_ptr   =
      SWIG_TypeQuery("NOX::Solver::TensorBased*");
    // Build a NOX::Solver::Generic object via the buildSolver factory
    Teuchos::RCP<NOX::Solver::Generic> rcp_solver =
      NOX::Solver::buildSolver(grp, tests, params);
    Teuchos::Ptr<NOX::Solver::Generic> solver_ptr = rcp_solver.release();
    // Try to downcast to a derived class
    {
      NOX::Solver::LineSearchBased * result =
	reinterpret_cast<NOX::Solver::LineSearchBased*>(solver_ptr.get());
      if (result)
	return SWIG_NewPointerObj(result, swig_NSLSB_ptr, 1);
    }
    {
      NOX::Solver::TrustRegionBased * result =
	reinterpret_cast<NOX::Solver::TrustRegionBased*>(solver_ptr.get());
      if (result)
	return SWIG_NewPointerObj(result, swig_NSTRB_ptr, 1);
    }
    {
      NOX::Solver::InexactTrustRegionBased * result =
	reinterpret_cast<NOX::Solver::InexactTrustRegionBased*>(solver_ptr.get());
      if (result)
	return SWIG_NewPointerObj(result, swig_NSITRB_ptr, 1);
    }
    {
      NOX::Solver::TensorBased * result =
	reinterpret_cast<NOX::Solver::TensorBased*>(solver_ptr.get());
      if (result)
	return SWIG_NewPointerObj(result, swig_NSTB_ptr, 1);
    }
    PyErr_SetString(PyExc_RuntimeError, "NOX::Solver::buildSolver returned unrecognized "
		    "derivative of NOX::Solver::Generic");
    return NULL;
  }
}


// NOX default solver.  Provide a function that returns a solver with
// default settings.
// %pythoncode
// {
// def defaultSolver(group, myPID=0, maxIt=100, tol=1.0e-4, precMaxAge=5):

//     # Convergence status test
//     nlParams = {"Nonlinear Solver" : "Line Search Based",
//                 "Printing"         : {"MyPID"            : myPID,
//                                       "Output Precision" : 3,
//                                       "Output Processor" : 0    },
//                 "Line Search"      : {"Method" : "Full Step"},
//                 "Direction"        : {"Method" : "Newton"},
//                 "Newton"           : {"Forcing Term Method" : "Constant"},
//                 "Linear Solver"    : {"Aztec Solver"    : "GMRES",
//                                       "Max Iterations"  : maxIt,
//                                       "Tolerance"       : tol,
//                                       "Preconditioner"  : "Ifpack",
//                                       "Max Age Of Prec" : precMaxAge    },
//                 "Solver Options"   : {"Status Test Check Type" : "Complete"}
//                 }
// }
    
// Turn off the exception handling
%exception;
