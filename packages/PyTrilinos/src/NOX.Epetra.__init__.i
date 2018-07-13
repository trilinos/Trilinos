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

%define %nox_epetra_docstring
"
PyTrilinos.NOX.Epetra is the python interface to namespace Epetra for
the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Epetra is to provide a concrete interface beteen
NOX and Epetra.

NOX.Epetra provides the following user-level classes:

    * Group                    - Epetra implementation of Abstract.Group
    * Vector                   - Epetra implementation of Abstract.Vector
    * FiniteDifference         - Class for estimating Jacobian w/finite differences
    * FiniteDifferenceColoring - FiniteDifference class, w/coloring efficiencies
    * MatrixFree               - Base class for Jacobian-free algorithms
    * Scaling                  - Class for controlling scalling of algebraic objects
    * LinearSystem             - Base class for interface to linear solvers
    * LinearSystemAztecOO      - Concrete implementation of LinearSystem
"
%enddef

%define %nox_epetra_import_code
"
from . import ___init__
"
%enddef

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
        moduleimport = %nox_epetra_import_code,
	docstring    = %nox_epetra_docstring) __init__

%{
// System include files
#include <vector>

// Configuration
#include "PyTrilinos_config.h"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#include "PyTrilinos_Epetra_Headers.hpp"

// EpetraExt include files
#ifdef HAVE_NOX_EPETRAEXT
#include "PyTrilinos_EpetraExt_Headers.hpp"
#endif

// NOX include files
#include "PyTrilinos_NOX_Abstract_Headers.hpp"
#include "PyTrilinos_NOX_Epetra_Headers.hpp"

// NumPy include
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Namespace flattening
using Teuchos::RCP;
using Teuchos::rcp;
using namespace NOX;
using namespace NOX::Abstract;
using namespace NOX::Epetra;
%}

// Configuration
%include "Epetra_DLLExportMacro.h"

// Include NOX documentation
%include "NOX_dox.i"

// General ignore directives
%ignore *::print(ostream &);
%ignore *::print(std::ostream &) const;
%ignore *::print(std::ostream &, int) const;
%ignore *::operator=;
%ignore *::operator<<;
%ignore *::operator[];

// SWIG library include files
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"
%teuchos_rcp(NOX::Abstract::Group)
%teuchos_rcp(NOX::Epetra::Interface::Required)
%teuchos_rcp(NOX::Epetra::Interface::Jacobian)
%teuchos_rcp(NOX::Epetra::Interface::Preconditioner)

// Allow import from the this directory and its parent, and force
// correct import of ___init__
%pythoncode
%{
import sys, os.path as op
thisDir   = op.dirname(op.abspath(__file__))
parentDir = op.normpath(op.join(thisDir,".."))
if not thisDir   in sys.path: sys.path.append(thisDir)
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}

// Include typemaps for Abstract base classes
%ignore *::getXPtr;
%ignore *::getFPtr;
%ignore *::getGradientPtr;
%ignore *::getNewtonPtr;
%include "NOX.Abstract_typemaps.i"
%import(module="Abstract" ) "NOX_Abstract_Group.H"
%import(module="Abstract" ) "NOX_Abstract_PrePostOperator.H"
%import(module="Abstract" ) "NOX_Abstract_MultiVector.H"
%import(module="Abstract" ) "NOX_Abstract_Vector.H"
%import(module="Interface") "NOX_Epetra_Interface_Required.H"
%import(module="Interface") "NOX_Epetra_Interface_Jacobian.H"
%import(module="Interface") "NOX_Epetra_Interface_Preconditioner.H"

// Support for Teuchos::RCPs
%teuchos_rcp(NOX::Epetra::Group)
%teuchos_rcp(NOX::Epetra::FiniteDifference)
%teuchos_rcp(NOX::Epetra::MatrixFree)
%teuchos_rcp(NOX::Epetra::LinearSystem)
%teuchos_rcp(NOX::Epetra::LinearSystemAztecOO)
%teuchos_rcp(NOX::Epetra::Scaling)
%teuchos_rcp(NOX::Epetra::VectorSpace)

// Epetra import
%import "Epetra.i"
%#if PY_VERSION_HEX >= 0x03000000

// EpetraExt import
#ifdef HAVE_NOX_EPETRAEXT
%ignore EpetraExt::Add;
%import "EpetraExt.i"
#endif

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

//////////////////////////////
// NOX.Epetra.Group support //
//////////////////////////////
%rename(Group_None) NOX::Epetra::Group::None;
%include "NOX_Epetra_Group.H"

///////////////////////////////
// NOX.Epetra.Vector support //
///////////////////////////////
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);
%ignore NOX::Epetra::Vector::getEpetraVector() const;
%include "NOX_Epetra_Vector.H"

/////////////////////////////////////////
// NOX.Epetra.FiniteDifference support //
/////////////////////////////////////////
%include "NOX_Epetra_FiniteDifference.H"

/////////////////////////////////////////////////
// NOX.Epetra.FiniteDifferenceColoring support //
/////////////////////////////////////////////////
#ifdef HAVE_NOX_EPETRAEXT
%teuchos_rcp(NOX::Epetra::FiniteDifferenceColoring)
namespace NOX
{
namespace Epetra
{
%extend FiniteDifferenceColoring
{
  FiniteDifferenceColoring(Teuchos::ParameterList & printingParams,
			   const Teuchos::RCP< Interface::Required > & i,
			   const NOX::Epetra::Vector & initialGuess,
			   const Teuchos::RCP< Epetra_CrsGraph > & rawGraph,
			   bool parallelColoring = false,
			   bool distance1 = false,
			   double beta = 1.0e-6,
			   double alpha = 1.0e-4)
  {
    // Construct the coloring algorithm functor and generate the color map
    EpetraExt::CrsGraph_MapColoring *mapColor = new EpetraExt::CrsGraph_MapColoring();
    const Teuchos::RCP< Epetra_MapColoring > colorMap = Teuchos::rcp(&(*mapColor)(*rawGraph));
    // Construct the color index functor and generate the column indexes
    EpetraExt::CrsGraph_MapColoringIndex *mapColorIndex =
      new EpetraExt::CrsGraph_MapColoringIndex(*colorMap);
    const Teuchos::RCP< std::vector< Epetra_IntVector > > columns =
      Teuchos::rcp(&(*mapColorIndex)(*rawGraph));
    // Construct the FiniteDifferenceColoring object
    FiniteDifferenceColoring *fdc =
      new FiniteDifferenceColoring(printingParams, i, initialGuess, rawGraph, colorMap,
				   columns, parallelColoring, distance1, beta, alpha);
    // Delete temporary functors
    delete mapColor;
    delete mapColorIndex;
    // Return the pointer to FiniteDifferenceColoring object
    return fdc;
  }
}
%ignore FiniteDifferenceColoring::FiniteDifferenceColoring;
}
}
%include "NOX_Epetra_FiniteDifferenceColoring.H"
#endif

///////////////////////////////////
// NOX.Epetra.MatrixFree support //
///////////////////////////////////
%include "NOX_Epetra_MatrixFree.H"

////////////////////////////////
// NOX.Epetra.Scaling support //
////////////////////////////////
%ignore operator<<(ostream&, NOX::Epetra::Scaling&);
%rename(Scaling_None) NOX::Epetra::Scaling::None;
%include "NOX_Epetra_Scaling.H"

/////////////////////////////////////
// NOX.Epetra.LinearSystem support //
/////////////////////////////////////
%feature("director") NOX::Epetra::LinearSystem;
// The following #define is to change the name of NOX method
// arguments that conflict with a SWIG director method argument
#define result nox_result
%include "NOX_Epetra_LinearSystem.H"
#undef result

////////////////////////////////////////////
// NOX.Epetra.LinearSystemAztecOO support //
////////////////////////////////////////////
// There are two LinearSystemAztecOO constructors that are similar,
// but different.  SWIG gets confused and thinks one overshadows the
// other.  This suppresses the warning.
%warnfilter(509) NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO;
%include "NOX_Epetra_LinearSystem_AztecOO.H"

////////////////////////////////////////////////
// NOX.Epetra.ModelEvaluatorInterface support //
////////////////////////////////////////////////
#ifdef HAVE_NOX_EPETRAEXT
%teuchos_rcp(NOX::Epetra::ModelEvaluatorInterface)
%import "EpetraExt.i"
namespace NOX
{
namespace Epetra
{
%extend ModelEvaluatorInterface
{
  ModelEvaluatorInterface(EpetraExt::ModelEvaluator & eeme)
  {
    Teuchos::RCP<EpetraExt::ModelEvaluator> eeme_ptr = Teuchos::rcpFromRef(eeme);
    return new ModelEvaluatorInterface(eeme_ptr);
  }
}
// Only accept the above constructor for wrapping
%ignore ModelEvaluatorInterface::ModelEvaluatorInterface;
}
}
%include "NOX_Epetra_ModelEvaluatorInterface.H"
#endif

// Turn off the exception handling
%exception;

///////////////////////
// Default factories //
///////////////////////

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
        comm = PyTrilinos.Epetra.SerialComm()
    myPID = comm.MyPID()

    # Create the printing parameter list
    Utils = PyTrilinos.NOX.Utils
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
    reqInterface        - a NOX.Epetra.Interface.Required object that defines
                          the interface to the nonlinear problem.  May be
                          None if both the Jacobian and preconditioner are
                          provided.
    jacInterface        - a NOX.Epetra.Interface.Jacobian object that defines the
                          Jacobian of the nonlinear problem.  Default None.
    jacobian            - if jacInterface is provided, this is the Epetra.Operator
                          that defines the Jacobian matrix.  Default None.
    precInterface       - a NOX.Epetra.Interface.Preconditioner object that defines
                          the preconditioner to the nonlinear problem.  Default None.
    preconditioner      - if precInterface is provided, this is the
                          Epetra.Operator that defines the preconditioner.
                          Default None.
    """

    # Sanity checks to prevent more cryptic problems down the road...
    if jacInterface is None or precInterface is None:
        assert isinstance(reqInterface, Interface.Required)
    if jacInterface is not None:
        assert isinstance(jacInterface, Interface.Jacobian        )
        assert isinstance(jacobian    , PyTrilinos.Epetra.Operator)
    if precInterface is not None:
        assert isinstance(precInterface , Interface.Preconditioner  )
        assert isinstance(preconditioner, PyTrilinos.Epetra.Operator)

    # Extract parameter lists
    printParams = nonlinearParameters["Printing"     ]
    lsParams    = nonlinearParameters["Linear Solver"]

    # Construct a NOX.Epetra.Vector from the Epetra.Vector
    clone = Vector(initGuess, Vector.CreateView)

    # Construct the linear system
    if jacInterface:
        if precInterface:
            linSys = LinearSystemAztecOO(printParams, lsParams,
                                         jacInterface, jacobian,
                                         precInterface, preconditioner,
                                         clone)
        else:
            linSys = LinearSystemAztecOO(printParams, lsParams,
                                         reqInterface,
                                         jacInterface, jacobian,
                                         clone)
    else:
        if precInterface:
            linSys = LinearSystemAztecOO(printParams, lsParams,
                                         reqInterface,
                                         precInterface, preconditioner,
                                         clone)
        else:
            linSys = LinearSystemAztecOO(printParams, lsParams,
                                         reqInterface,
                                         clone)

    # Construct and return the default Group
    group = Group(printParams, reqInterface, clone, linSys)
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
    StatusTest = PyTrilinos.NOX.StatusTest
    converged  = StatusTest.Combo(StatusTest.Combo.AND)
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
                  precInterface=None, preconditioner=None, nlParams=None,
                  absTol=1.0e-8, relTol=1.0e-2, relGroup=None, updateTol=1.0e-5,
                  wAbsTol=1.0e-8, wRelTol=1.0e-2, maxIters=20, finiteValue=True):
    """
    defaultSolver(initGuess, reqInterface, jacInterface=None, jacobian=None,
                  precInterface=None, preconditioner=None, nlParams=None) -> Solver

    Return a default NOX Solver based on the given arguments:

    initGuess      - an Epetra.Vector initial guess
    reqInterface   - a NOX.Epetra.Interface.Required object that defines
                     the interface to the nonlinear problem.  May be
                     None if both the Jacobian and preconditioner are
                     provided.
    jacInterface   - a NOX.Epetra.Interface.Jacobian object that defines the
                     Jacobian of the nonlinear problem.  Default None.
    jacobian       - if jacInterface is provided, this should also be provided
                     and is the Epetra.Operator that defines the Jacobian
                     matrix.  Default None. 
    precInterface  - a NOX.Epetra.Interface.Preconditioner object that defines
                     the preconditioner to the nonlinear problem.  Default None.
    preconditioner - if precInterface is provided, this should also be provided
                     and is the Epetra.Operator that defines the preconditioner.
                     Default None.
    nlParams       - dict that contains a list of nonlinear parameters.  Default
                     None, in which case defaultNonlinearParameters() is used.
    absTol         - if not None, include an absolute residual status test,
                     using this value as the tolerance.  Default 1.0e-8.
    relTol         - if not None, include a relative residual status test, using
                     this value as the tolerance.  Default 1.0e-2.
    relGroup       - if relTol is specified, use this Group to determine the
                     scaling.  If relGroup is None, use the result of
                     defaultGroup().  Default None.
    updateTol      - if not None, include an update status test, using this
                     value as the tolerance.  Default 1.0e-5.
    wAbsTol        - if not None, along with wRelTol, include a weighted RMS
                     status test, using this value as the absolute tolerance.
                     Default 1.0e-8.
    wRelTol        - if not None, along with wAbsTol, include a weighted RMS
                     status test, using this value as the relative tolerance.
                     Default 1.0e-2.
    maxIters       - if not None, include a maximum nonlinear iterations status
                     test, using this value as the maximum allowable iterations.
                     Default 20. 
    finiteValue    - if True, include a finite value status test.  Default
                     True. 
    """

    # Sanity checks to prevent more cryptic problems down the road...
    if jacInterface is None or precInterface is None:
        assert isinstance(reqInterface, Interface.Required)
    if jacInterface is not None:
        assert isinstance(jacInterface, Interface.Jacobian        )
        assert isinstance(jacobian    , PyTrilinos.Epetra.Operator)
    if precInterface is not None:
        assert isinstance(precInterface , Interface.Preconditioner  )
        assert isinstance(preconditioner, PyTrilinos.Epetra.Operator)

    # Get the communicator
    comm = initGuess.Comm()

    # Get the nonlinear parameters
    if nlParams is None:
        nlParams = defaultNonlinearParameters(comm,2)

    # Build the default Group
    group = defaultGroup(nlParams, initGuess, reqInterface, jacInterface,
                         jacobian, precInterface, preconditioner)

    # Get the default StatusTest
    if relTol and (relGroup is None): relGroup = group
    statusTest = defaultStatusTest(absTol      = absTol,
                                   relTol      = relTol,
                                   relGroup    = relGroup,
                                   updateTol   = updateTol,
                                   wAbsTol     = wAbsTol,
                                   wRelTol     = wRelTol,
                                   maxIters    = maxIters,
                                   finiteValue = finiteValue)

    # Return the default Solver
    solver = PyTrilinos.NOX.Solver.buildSolver(group, statusTest, nlParams)
    #solver.group      = group        ### By adding group, statusTest and
    solver.statusTest = statusTest   ### nlParams as attributes to the Solver
    solver.nlParams   = nlParams     ### variable, we ensure that they do not
                                     ### get destroyed.  This is a workaround for
                                     ### a Teuchos::RCP wrapper bug.
    return solver
%}
