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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %nox_petsc_docstring
"
PyTrilinos.NOX.PETSc is the python interface to namespace PETSc for
the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.PETSc is to provide a concrete interface beteen
NOX and PETSc.

NOX.PETSC provides the following user-level classes:

    * Group                    - PETSc implementation of Abstract.Group
    * Vector                   - PETSc implementation of Abstract.Vector
    * SharedJacobian           - PETSc implementation of SharedJacobian
"
%enddef

%module(package      = "PyTrilinos.NOX.PETSc",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_petsc_docstring) __init__

%{
// Configuration
#include "PyTrilinos_config.h"

// Teuchos includes
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

// NOX includes
#include "NOX_Petsc.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Namespace flattening
using namespace NOX::Petsc;
%}

// Include the PETSc4Py SWIG interface file
%include "petsc4py/petsc4py.i"

// General exception handling
%include "exception.i"

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

// Include NOX documentation
%include "NOX_dox.i"

// General ignore directives
%ignore *::print(ostream &);
%ignore *::print(std::ostream &) const;
%ignore *::print(std::ostream &, int) const;
%ignore *::operator=;
%ignore *::operator<<;
%ignore *::operator[];

// SWIG library includes
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"

// Support for Teuchos::RCPs

// *** Don't erase !!! ***
// Include typemaps for converting raw types to NOX.Abstract types
//%include "NOX.Abstract_typemaps.i"
// ***

// Allow import from the parent directory
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}

// NOX base classes
%teuchos_rcp(NOX::Abstract::Group)
%import(module="Abstract") "NOX_Abstract_Group.H"
%import(module="Abstract") "NOX_Abstract_PrePostOperator.H"
%import(module="Abstract") "NOX_Abstract_MultiVector.H"
%import(module="Abstract") "NOX_Abstract_Vector.H"

// NOX::Petsc::Interface imports
%teuchos_rcp(NOX::Petsc::Interface)
%import(module="Interface") "NOX_Petsc_Interface.H"

/////////////////////////////
// NOX.Petsc.Group support //
/////////////////////////////
%teuchos_rcp(NOX::Petsc::Group)
%rename(Group_None) NOX::Petsc::Group::None;
%include "NOX_Petsc_Group.H"

//////////////////////////////
// NOX.Petsc.Vector support //
//////////////////////////////
%ignore NOX::Petsc::Vector::Vector(const Vec &, std::string,  NOX::CopyType);
%include "NOX_Petsc_Vector.H"

//////////////////////////////////////
// NOX.Petsc.SharedJacobian support //
//////////////////////////////////////
%teuchos_rcp(NOX::Petsc::SharedJacobian)
%include "NOX_Petsc_SharedJacobian.H"

// Turn off the exception handling
%exception;

// ///////////////////////
// // Default factories //
// ///////////////////////

// // defaultSolver() and supporting functions
// %pythoncode
// %{
// def defaultNonlinearParameters(comm=None, verbosity=0, outputPrec=3,
//                                maxIterations=800, tolerance=1.0e-4):
//     """
//     defaultNonlinearParameters(comm=None, verbosity=0, outputPrec=3,
//                                maxIterations=800, tolerance=1.0e-4) -> dict

//     Return a dictionary that can serve as a default list of parameters for a NOX
//     solver.  Entries can be altered before passing to a NOX solver constructor.

//     comm          - Epetra communicator object.  If not provided, the function
//                     uses an Epetra.SerialComm communicator.

//     verbosity     - A simple indication of verbosity level.  0: errors and test
//                     details.  1: debugging, warnings, details, parameters and
//                     linear solver details.  2: inner iteration, outer iteration
//                     status test and outer iteration.  Default 0.

//     outputPrec    - Number of significant digits to output.  Default 3.

//     maxIterations - Maximum allowable linear (inner) iterations.  Default 800.

//     tolerance     - Linear solver tolerance level.  Default 1.0e-4.
//     """
//     # Communicator
//     if comm is None:
//         comm = PyTrilinos.Epetra.SerialComm()
//     myPID = comm.MyPID()

//     # Create the printing parameter list
//     Utils = PyTrilinos.NOX.Utils
//     outputInfo = Utils.Error + Utils.TestDetails
//     if verbosity: outputInfo += Utils.Debug      + \
//                                 Utils.Warning    + \
//                                 Utils.Details    + \
//                                 Utils.Parameters + \
//                                 Utils.LinearSolverDetails
//     if verbosity > 1: outputInfo += Utils.InnerIteration           + \
//                                     Utils.OuterIterationStatusTest + \
//                                     Utils.OuterIteration
//     printParams = {"MyPID"              : myPID,
//                    "Output Precision"   : outputPrec,
//                    "Output Processor"   : 0,
//                    "Output Information" : outputInfo}

//     # Create the linear solver parameter list
//     lsParams = {"Aztec Solver"    : "GMRES",
//                 "Max Iterations"  : maxIterations,
//                 "Tolerance"       : tolerance,
//                 "Preconditioner"  : "Ifpack",
//                 "Max Age Of Prec" : 5       }

//     # Create the nonlinear solver parameter list
//     nlParams = {"Nonlinear Solver" : "Line Search Based",
//                 "Printing"         : printParams,
//                 "Line Search"      : {"Method" : "Full Step"},
//                 "Direction"        : {"Method" : "Newton"},
//                 "Newton"           : {"Forcing Term Method" : "Constant"},
//                 "Linear Solver"    : lsParams,
//                 "Solver Options"   : {"Status Test Check Type" : "Complete"}
//                 }
    
//     return nlParams

// def defaultGroup(nonlinearParameters, initGuess, reqInterface, jacInterface=None,
//                  jacobian=None, precInterface=None, preconditioner=None):
//     """
//     defaultGroup(nonlinearParameters, initGuess, reqInterface, jacInterface=None,
//                  jacobian=None, precInterface=None, preconditioner=None) -> Group

//     Return a NOX.Epetra.Group based upon the given input arguments:

//     nonlinearParameters - a dict with nonlinear parameters.  Can be obtained
//                           from defaultNonlinearParameters()
//     initGuess           - an initial guess Epetra.Vector.
//     reqInterface        - a NOX.Epetra.Interface.Required object that defines
//                           the interface to the nonlinear problem.  May be
//                           None if both the Jacobian and preconditioner are
//                           provided.
//     jacInterface        - a NOX.Epetra.Interface.Jacobian object that defines the
//                           Jacobian of the nonlinear problem.  Default None.
//     jacobian            - if jacInterface is provided, this is the Epetra.Operator
//                           that defines the Jacobian matrix.  Default None.
//     precInterface       - a NOX.Epetra.Interface.Preconditioner object that defines
//                           the preconditioner to the nonlinear problem.  Default None.
//     preconditioner      - if precInterface is provided, this is the
//                           Epetra.Operator that defines the preconditioner.
//                           Default None.
//     """

//     # Sanity checks to prevent more cryptic problems down the road...
//     if jacInterface is None or precInterface is None:
//         assert isinstance(reqInterface, Interface.Required)
//     if jacInterface is not None:
//         assert isinstance(jacInterface, Interface.Jacobian        )
//         assert isinstance(jacobian    , (PyTrilinos.Epetra.Operator, Epetra.Operator))
//     if precInterface is not None:
//         assert isinstance(precInterface , Interface.Preconditioner  )
//         assert isinstance(preconditioner, (PyTrilinos.Epetra.Operator, Epetra.Operator))

//     # Extract parameter lists
//     printParams = nonlinearParameters["Printing"     ]
//     lsParams    = nonlinearParameters["Linear Solver"]

//     # Construct a NOX.Epetra.Vector from the Epetra.Vector
//     clone = Vector(initGuess, Vector.CreateView)

//     # Construct the linear system
//     if jacInterface:
//         if precInterface:
//             linSys = LinearSystemAztecOO(printParams, lsParams,
//                                          jacInterface, jacobian,
//                                          precInterface, preconditioner,
//                                          clone)
//         else:
//             linSys = LinearSystemAztecOO(printParams, lsParams,
//                                          reqInterface,
//                                          jacInterface, jacobian,
//                                          clone)
//     else:
//         if precInterface:
//             linSys = LinearSystemAztecOO(printParams, lsParams,
//                                          reqInterface,
//                                          precInterface, preconditioner,
//                                          clone)
//         else:
//             linSys = LinearSystemAztecOO(printParams, lsParams,
//                                          reqInterface,
//                                          clone)

//     # Construct and return the default Group
//     group = Group(printParams, reqInterface, clone, linSys)
//     group.linSys = linSys   ### By adding linSys as an attribute to the Group
//                             ### variable, we ensure that linSys does not get
//                             ### destroyed.  This is a workaround for a
//                             ### Teuchos::RCP wrapper bug.
//     return group

// def defaultStatusTest(absTol=None, relTol=None, relGroup=None, updateTol=None,
//                       wAbsTol=None, wRelTol=None, maxIters=None,
//                       finiteValue=False):
//     """
//     defaultStatusTest(absTol=None, relTol=None, relGroup=None, updateTol=None,
//                       wAbsTol=None, wRelTol=None, maxIters=None,
//                       finiteValue=False) -> StatusTest

//     Return a StatusTest object based upon the input arguments:

//     absTol      - if specified, include an absolute residual status test, using
//                   this value as the tolerance
//     relTol      - if specified, along with relGroup, include a relative residual
//                   status test, using this value as the tolerance
//     relGroup    - if specified, along with relTol, include a relative residual
//                   status test, using this Group to determine the scaling
//     updateTol   - if specified, include an update status test, using this value
//                   as the tolerance
//     wAbsTol     - if specified, along with wRelTol, include a weighted RMS
//                   status test, using this value as the absolute tolerance
//     wRelTol     - if specified, along with wAbsTol, include a weighted RMS
//                   status test, using this value as the relative tolerance
//     maxIters    - if specified, include a maximum iterations status test, using
//                   this value as the maximum allowable iterations
//     finiteValue - if True, include a finite value status test.  Default False.
//     """
//     # Build the convergence portion of the status test
//     StatusTest = PyTrilinos.NOX.StatusTest
//     converged  = StatusTest.Combo(StatusTest.Combo.AND)
//     converged.tests = [ ]   ### By adding this list of tests as an attribute to
//                             ### the StatusTest variables, we ensure that linSys
//                             ### does not get destroyed.  This is a workaround
//                             ### for a Teuchos::RCP wrapper bug.

//     if absTol:
//         absTest = StatusTest.NormF(absTol)
//         converged.addStatusTest(absTest)
//         converged.tests.append(absTest)
//     if relGroup and relTol:
//         relTest = StatusTest.NormF(relGroup,relTol)
//         converged.addStatusTest(relTest)
//         converged.tests.append(relTest)
//     if wAbsTol and wRelTol:
//         wrmsTest = StatusTest.NormWRMS(wRelTol,wAbsTol)
//         converged.addStatusTest(wrmsTest)
//         converged.tests.append(wrmsTest)
//     if updateTol:
//         updateTest = StatusTest.NormUpdate(updateTol)
//         converged.addStatusTest(updateTest)
//         converged.tests.append(updateTest)

//     # Return if done
//     if not (maxIters or finiteValue):
//         return converged

//     # Add the divergence portion of the default status test
//     combo = StatusTest.Combo(StatusTest.Combo.OR)
//     combo.tests = [ ]
//     if finiteValue:
//         fvTest = StatusTest.FiniteValue()
//         combo.addStatusTest(fvTest)
//         combo.tests.append(fvTest)
//     combo.addStatusTest(converged)
//     combo.tests.append(converged)
//     if maxIters:
//         maxIterTest = StatusTest.MaxIters(maxIters)
//         combo.addStatusTest(maxIterTest)
//         combo.tests.append(maxIterTest)
//     return combo

// def defaultSolver(initGuess, reqInterface, jacInterface=None, jacobian=None,
//                   precInterface=None, preconditioner=None, nlParams=None,
//                   absTol=1.0e-8, relTol=1.0e-2, relGroup=None, updateTol=1.0e-5,
//                   wAbsTol=1.0e-8, wRelTol=1.0e-2, maxIters=20, finiteValue=True):
//     """
//     defaultSolver(initGuess, reqInterface, jacInterface=None, jacobian=None,
//                   precInterface=None, preconditioner=None, nlParams=None) -> Solver

//     Return a default NOX Solver based on the given arguments:

//     initGuess      - an Epetra.Vector initial guess
//     reqInterface   - a NOX.Epetra.Interface.Required object that defines
//                      the interface to the nonlinear problem.  May be
//                      None if both the Jacobian and preconditioner are
//                      provided.
//     jacInterface   - a NOX.Epetra.Interface.Jacobian object that defines the
//                      Jacobian of the nonlinear problem.  Default None.
//     jacobian       - if jacInterface is provided, this should also be provided
//                      and is the Epetra.Operator that defines the Jacobian
//                      matrix.  Default None. 
//     precInterface  - a NOX.Epetra.Interface.Preconditioner object that defines
//                      the preconditioner to the nonlinear problem.  Default None.
//     preconditioner - if precInterface is provided, this should also be provided
//                      and is the Epetra.Operator that defines the preconditioner.
//                      Default None.
//     nlParams       - dict that contains a list of nonlinear parameters.  Default
//                      None, in which case defaultNonlinearParameters() is used.
//     absTol         - if not None, include an absolute residual status test,
//                      using this value as the tolerance.  Default 1.0e-8.
//     relTol         - if not None, include a relative residual status test, using
//                      this value as the tolerance.  Default 1.0e-2.
//     relGroup       - if relTol is specified, use this Group to determine the
//                      scaling.  If relGroup is None, use the result of
//                      defaultGroup().  Default None.
//     updateTol      - if not None, include an update status test, using this
//                      value as the tolerance.  Default 1.0e-5.
//     wAbsTol        - if not None, along with wRelTol, include a weighted RMS
//                      status test, using this value as the absolute tolerance.
//                      Default 1.0e-8.
//     wRelTol        - if not None, along with wAbsTol, include a weighted RMS
//                      status test, using this value as the relative tolerance.
//                      Default 1.0e-2.
//     maxIters       - if not None, include a maximum nonlinear iterations status
//                      test, using this value as the maximum allowable iterations.
//                      Default 20. 
//     finiteValue    - if True, include a finite value status test.  Default
//                      True. 
//     """

//     # Sanity checks to prevent more cryptic problems down the road...
//     if jacInterface is None or precInterface is None:
//         assert isinstance(reqInterface, Interface.Required)
//     if jacInterface is not None:
//         assert isinstance(jacInterface, Interface.Jacobian        )
//         assert isinstance(jacobian    , (PyTrilinos.Epetra.Operator, Epetra.Operator))
//     if precInterface is not None:
//         assert isinstance(precInterface , Interface.Preconditioner  )
//         assert isinstance(preconditioner, (PyTrilinos.Epetra.Operator, Epetra.Operator))

//     # Get the communicator
//     comm = initGuess.Comm()

//     # Get the nonlinear parameters
//     if nlParams is None:
//         nlParams = defaultNonlinearParameters(comm,2)

//     # Build the default Group
//     group = defaultGroup(nlParams, initGuess, reqInterface, jacInterface,
//                          jacobian, precInterface, preconditioner)

//     # Get the default StatusTest
//     if relTol and (relGroup is None): relGroup = group
//     statusTest = defaultStatusTest(absTol      = absTol,
//                                    relTol      = relTol,
//                                    relGroup    = relGroup,
//                                    updateTol   = updateTol,
//                                    wAbsTol     = wAbsTol,
//                                    wRelTol     = wRelTol,
//                                    maxIters    = maxIters,
//                                    finiteValue = finiteValue)

//     # Return the default Solver
//     solver = PyTrilinos.NOX.Solver.buildSolver(group, statusTest, nlParams)
//     #solver.group      = group        ### By adding group, statusTest and
//     solver.statusTest = statusTest   ### nlParams as attributes to the Solver
//     solver.nlParams   = nlParams     ### variable, we ensure that they do not
//                                      ### get destroyed.  This is a workaround for
//                                      ### a Teuchos::RCP wrapper bug.
//     return solver
// %}
