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

%define %nox_docstring
"
PyTrilinos.NOX is the python interface to the Trilinos nonlinear
solver package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX is to provide robust nonlinear solvers for the
problem of finding x such that F(x)=0.  In C++, NOX supports several
namespaces, some of which are sub-modules in python:

    * Abstract          - Base classes for abstract interface to NOX
    * Epetra            - Epetra implementation
    * Epetra.Interface  - Concrete interface for Epetra
    * Solver            - Solver manager class and supporting utilities
    * StatusTest        - Support for customizable stopping criteria

The top-level NOX module provides the following user-level class:

    * Utils  - Various utilities

For an example of usage of all of NOX, please consult the following
script in the example subdirectory of the PyTrilinos package:

    * exNOX_1Dfem.py
"
%enddef

%module(package   = "PyTrilinos.NOX",
	autodoc   = "1",
	docstring = %nox_docstring) __init__

%{
// System includes
#include <sstream>

// PyTrilinos configuration
#include "PyTrilinos_config.h"

// Teuchos include
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_Version.H"
#include "NOX_Utils.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// General ignore directives
%ignore operator<<;
%ignore *::operator=;

// Auto-documentation feature
%feature("autodoc", "1");

// Include NOX documentation
%include "NOX_dox.i"

// SWIG library includes
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"
// Note: Teuchos.i turns off warnings for nested classes, so we do not
// have to do it again.

//////////////////////////////////////
// PyTrilinos configuration support //
//////////////////////////////////////
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%constant bool Have_Epetra = true;
#else
%constant bool Have_Epetra = false;
#endif

/////////////////////////
// NOX Version support //
/////////////////////////
%include "NOX_Version.H"
%pythoncode
%{
__version__ = version().split()[2]
%}

///////////////////////
// NOX Utils support //
///////////////////////
%rename(_print) NOX::Utils::print;
%ignore NOX::Utils::fill;
%ignore NOX::Utils::sciformat;
%include "NOX_Utils.H"

// NOX namespace imports
%pythoncode
%{
# Abstract, Solver, and StatusTest namespaces
__all__ = ['Abstract', 'Solver', 'StatusTest']
import Abstract
import Solver
import StatusTest
%}

// NOX.Epetra namespace
#ifdef HAVE_NOX_EPETRA
%pythoncode
%{

# Epetra namespace
__all__.append('Epetra')
%}
#endif

// defaultSolver() and supporting functions
%pythoncode
%{
def defaultNonlinearParameters(comm=None, verbosity=0, outputPrec=3,
                               maxIterations=800, tolerance=1.0e-4):
    """
    defaultNonlinearParameters(comm=None, verbosity=0, outputPrec=3,
                               maxIterations=800, tolerance=1.0e-4) -> dict

    Return a dictionary that can serve as a default list of parameters for a NOX
    solver.  Entries can be altered before passing to a NOX solver constructor.

    comm          - Epetra communicator object.  If not provided, the function
                    uses an Epetra.SerialComm communicator.

    verbosity     - A simple indication of verbosity level.  0: errors and test
                    details.  1: debugging, warnings, details, parameters and
                    linear solver details.  2: inner iteration, outer iteration
                    status test and outer iteration.  Default 0.

    outputPrec    - Number of significant digits to output.  Default 3.

    maxIterations - Maximum allowable linear (inner) iterations.  Default 800.

    tolerance     - Linear solver tolerance level.  Default 1.0e-4.
    """
    # Communicator
    if comm is None:
        comm = Epetra.SerialComm()
    myPID = comm.MyPID()

    # Create the printing parameter list
    outputInfo = Utils.Error + Utils.TestDetails
    if verbosity: outputInfo += Utils.Debug      + \
                                Utils.Warning    + \
                                Utils.Details    + \
                                Utils.Parameters + \
                                Utils.LinearSolverDetails
    if verbosity > 1: outputInfo += Utils.InnerIteration           + \
                                    Utils.OuterIterationStatusTest + \
                                    Utils.OuterIteration
    printParams = {"MyPID"              : myPID,
                   "Output Precision"   : outputPrec,
                   "Output Processor"   : 0,
                   "Output Information" : outputInfo}

    # Create the linear solver parameter list
    lsParams = {"Aztec Solver"    : "GMRES",
                "Max Iterations"  : maxIterations,
                "Tolerance"       : tolerance,
                "Preconditioner"  : "Ifpack",
                "Max Age Of Prec" : 5       }

    # Create the nonlinear solver parameter list
    nlParams = {"Nonlinear Solver" : "Line Search Based",
                "Printing"         : printParams,
                "Line Search"      : {"Method" : "Full Step"},
                "Direction"        : {"Method" : "Newton"},
                "Newton"           : {"Forcing Term Method" : "Constant"},
                "Linear Solver"    : lsParams,
                "Solver Options"   : {"Status Test Check Type" : "Complete"}
                }
    
    return nlParams

def defaultGroup(nonlinearParameters, initGuess, reqInterface, jacInterface=None,
                 jacobian=None, precInterface=None, preconditioner=None):
    """
    defaultGroup(nonlinearParameters, initGuess, reqInterface, jacInterface=None,
                 jacobian=None, precInterface=None, preconditioner=None) -> Group

    Return a NOX.Epetra.Group based upon the given input arguments:

    nonlinearParameters - a dict with nonlinear parameters.  Can be obtained
                          from defaultNonlinearParameters()
    initGuess           - an initial guess Epetra.Vector.
    reqInterface        - an Interface.Required object
    jacInterface        - an Interface.Jacobian object
    jacobian            - if jacInterface is given, this should be the
                          Epetra.Operator that is the Jacobian matrix
    precInterface       - an Interface.Preconditioner object
    preconditioner      - if precInterface is given, this should be the
                          Epetra.Operator that is the Preconditioner
    """
    # Extract parameter lists
    printParams = nonlinearParameters["Printing"     ]
    lsParams    = nonlinearParameters["Linear Solver"]

    # Construct a NOX.Epetra.Vector from the Epetra.Vector
    clone = Epetra.Vector(initGuess, Epetra.Vector.CreateView)

    # Construct the linear system
    if jacInterface:
        if precInterface:
            linSys = Epetra.LinearSystemAztecOO(printParams, lsParams,
                                                jacInterface, jacobian,
                                                precInterface, preconditioner,
                                                clone)
        else:
            linSys = Epetra.LinearSystemAztecOO(printParams, lsParams,
                                                reqInterface,
                                                jacInterface, jacobian,
                                                clone)
    else:
        if precInterface:
            linSys = Epetra.LinearSystemAztecOO(printParams, lsParams,
                                                reqInterface,
                                                precInterface, preconditioner,
                                                clone)
        else:
            linSys = Epetra.LinearSystemAztecOO(printParams, lsParams,
                                                reqInterface,
                                                clone)

    # Construct and return the default Group
    group = Epetra.Group(printParams, reqInterface, clone, linSys)
    group.linSys = linSys   ### By adding linSys as an attribute to the Group
                            ### variable, we ensure that linSys does not get
                            ### destroyed.  This is a workaround for a
                            ### Teuchos::RCP wrapper bug.
    return group

def defaultStatusTest(absTol=None, relTol=None, relGroup=None, updateTol=None,
                      wAbsTol=None, wRelTol=None, maxIters=None,
                      finiteValue=False):
    """
    defaultStatusTest(absTol=None, relTol=None, relGroup=None, updateTol=None,
                      wAbsTol=None, wRelTol=None, maxIters=None,
                      finiteValue=False) -> StatusTest

    Return a StatusTest object based upon the input arguments:

    absTol      - if specified, include an absolute residual status test, using
                  this value as the tolerance
    relTol      - if specified, along with relGroup, include a relative residual
                  status test, using this value as the tolerance
    relGroup    - if specified, along with relTol, include a relative residual
                  status test, using this Group to determine the scaling
    updateTol   - if specified, include an update status test, using this value
                  as the tolerance
    wAbsTol     - if specified, along with wRelTol, include a weighted RMS
                  status test, using this value as the absolute tolerance
    wRelTol     - if specified, along with wAbsTol, include a weighted RMS
                  status test, using this value as the relative tolerance
    maxIters    - if specified, include a maximum iterations status test, using
                  this value as the maximum allowable iterations
    finiteValue - if True, include a finite value status test.  Default False.
    """
    # Build the convergence portion of the status test
    converged = StatusTest.Combo(StatusTest.Combo.AND)
    converged.tests = [ ]   ### By adding this list of tests as an attribute to
                            ### the StatusTest variables, we ensure that linSys
                            ### does not get destroyed.  This is a workaround
                            ### for a Teuchos::RCP wrapper bug.

    if absTol:
        absTest = StatusTest.NormF(absTol)
        converged.addStatusTest(absTest)
        converged.tests.append(absTest)
    if relGroup and relTol:
        relTest = StatusTest.NormF(relGroup,relTol)
        converged.addStatusTest(relTest)
        converged.tests.append(relTest)
    if wAbsTol and wRelTol:
        wrmsTest = StatusTest.NormWRMS(wRelTol,wAbsTol)
        converged.addStatusTest(wrmsTest)
        converged.tests.append(wrmsTest)
    if updateTol:
        updateTest = StatusTest.NormUpdate(updateTol)
        converged.addStatusTest(updateTest)
        converged.tests.append(updateTest)

    # Return if done
    if not (maxIters or finiteValue):
        return converged

    # Add the divergence portion of the default status test
    combo = StatusTest.Combo(StatusTest.Combo.OR)
    combo.tests = [ ]
    if finiteValue:
        fvTest = StatusTest.FiniteValue()
        combo.addStatusTest(fvTest)
        combo.tests.append(fvTest)
    combo.addStatusTest(converged)
    combo.tests.append(converged)
    if maxIters:
        maxIterTest = StatusTest.MaxIters(maxIters)
        combo.addStatusTest(maxIterTest)
        combo.tests.append(maxIterTest)
    return combo

def defaultSolver(initGuess, reqInterface, jacInterface=None, jacobian=None,
                  precInterface=None, preconditioner=None, nlParams=None):
    """
    defaultSolver(initGuess, reqInterface, jacInterface=None, jacobian=None,
                  precInterface=None, preconditioner=None, nlParams=None) -> Solver

    Return a default NOX Solver based on the given arguments:

    initGuess      - an Epetra.Vector initial guess
    reqInterface   - a NOX.Epetra.Interface.Required object that defines the
                     interface to the nonlinear problem
    jacInterface   - a NOX.Epetra.Interface.Jacobian object that defines the
                     Jacobian of the nonlinear problem. Default None.
    jacobian       - an Epetra.Operator that defines the Jacobian matrix.
                     Default None.
    precInterface  - a NOX.Epetra.Interface.Preconditioner object that defines
                     the preconditioner to the nonlinear problem. Default None.
    preconditioner - an Epetra.Operator that defines the preconditioner.
                     Default None.
    nlParams       - dict that contains a list of nonlinear parameters.  Default
                     None, in which case defaultNonlinearParameters() is used.
    """
    # Get the communicator
    comm = initGuess.Comm()

    # Get the nonlinear parameters
    if nlParams is None:
        nlParams = defaultNonlinearParameters(comm,2)

    # Build the default Group
    group = defaultGroup(nlParams, initGuess, reqInterface, jacInterface,
                         jacobian, precInterface, preconditioner)

    # Get the default StatusTest
    statusTest = defaultStatusTest(absTol      = 1.0e-8,
                                   relTol      = 1.0e-2,
                                   relGroup    = group,
                                   updateTol   = 1.0e-5,
                                   wAbsTol     = 1.0e-8,
                                   wRelTol     = 1.0e-2,
                                   maxIters    = 20,
                                   finiteValue = True)

    # Return the default Solver
    solver = Solver.buildSolver(group, statusTest, nlParams)
    solver.group      = group        ### By adding group, statusTest and
    solver.statusTest = statusTest   ### nlParams as attributes to the Solver
    solver.nlParams   = nlParams     ### variable, we ensure that they do not
                                     ### get destroyed.  This is a workaround for
                                     ### a Teuchos::RCP wrapper bug.
    return solver
%}
