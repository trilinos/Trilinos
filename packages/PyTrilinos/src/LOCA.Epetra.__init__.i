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

%define %loca_epetra_docstring
"
PyTrilinos.LOCA.Epetra is the python interface to namespace Epetra of
the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra is to provide an extension of the
NOX.Epetra.Group to LOCA.  The python version of LOCA.Epetra supports
the following classes:

    * Group  - Extension of the NOX.Epetra.Group to LOCA
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
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// PyTrilinos includes
#include "PyTrilinos_Teuchos_Util.hpp"
#include "PyTrilinos_Epetra_Util.hpp"

// Local Epetra includes
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"

// Epetra includes
#include "Epetra_DLLExportMacro.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

// NOX includes
#include "NOX.H"
#include "NOX_Epetra_Group.H"

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
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

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

%include "Epetra_DLLExportMacro.h"

// Teuchos and Epetra interface support
%import "Teuchos.i"
%include "Epetra_Base.i"    // For PyExc_EpetraError
%import "Epetra.i"

// Teuchos RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::PhaseTransition::AbstractGroup)
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::BorderedSystem::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::Group)
%teuchos_rcp(LOCA::Homotopy::DeflatedGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::ExtendedGroup)
%teuchos_rcp(LOCA::Hopf::MinimallyAugmented::Constraint)
%teuchos_rcp(LOCA::Pitchfork::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Pitchfork::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::Abstract::Group)
%teuchos_rcp(LOCA::Abstract::TransposeSolveGroup)

// NOX interface support
%import "NOX.Abstract.i"
%import "NOX.Epetra.__init__.i"
%import "NOX.Epetra.Interface.i"

// Allow import from the parent directory
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
from .. import Abstract
%}

// LOCA base classes
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterface.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
%import(module="PhaseTransition") "LOCA_PhaseTransition_AbstractGroup.H"
%import(module="TimeDependent") "LOCA_TimeDependent_AbstractGroup.H"
%import(module="BorderedSystem") "LOCA_BorderedSystem_AbstractGroup.H"
%import(module="Homotopy") "LOCA_Homotopy_AbstractGroup.H"
%import(module="Homotopy") "LOCA_Homotopy_Group.H"
%import(module="Homotopy") "LOCA_Homotopy_DeflatedGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
%import(module="TurningPoint.MinimallyAugmented") "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%import(module="Hopf.MooreSpence") "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
%import(module="Hopf.MinimallyAugmented") "LOCA_Hopf_MinimallyAugmented_Constraint.H"
%import(module="Pitchfork.MooreSpence") "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%import(module="Pitchfork.MinimallyAugmented") "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
%import(module="Abstract") "LOCA_Abstract_Group.H"
%import(module="Abstract") "LOCA_Abstract_TransposeSolveGroup.H"

// The above %import(module="Abstract") ... directives can cause an
// "import Abstract" to appear in the .py file, causing Abstract to
// point to NOX.Abstract.  Force it back to LOCA.Abstract.  Also,
// ___init__ was pointing to Pitchfork/___init__.so (?!?), so I fix
// that, too.
%pythoncode
%{
del ___init__
from .. import Abstract
from .  import ___init__
%}

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

/////////////////////////
// LOCA Epetra support //
/////////////////////////

#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
%include "LOCA_Epetra.H"

//////////////////////////////
// LOCA.Epetra.Group support //
//////////////////////////////

// temporarily ignore conflict-causing constructor.  TODO: fix this issue
%ignore LOCA::Epetra::Group::Group(Teuchos::RCP< LOCA::GlobalData > const &,Teuchos::ParameterList &,Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree > const &,NOX::Epetra::Vector &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,LOCA::ParameterVector const &);

%teuchos_rcp(LOCA::Epetra::Group)
%include "LOCA_Epetra_Group.H"

%pythoncode
%{
def defaultContinuationParameters(comm=None,
                                  verbosity=0,
                                  outputPrec=3,
                                  maxIterations=800,
                                  tolerance=1.0e-4):
    """
    defaultContinuationParameters(comm=None,
                                  verbosity=0,
                                  outputPrec=3,
                                  maxIterations=800,
                                  tolerance=1.0e-4) -> dict

    Return a dictionary that can serve as a default list of parameters for LOCA
    constructors.  Entries can be altered before passing to a LOCA constructor.

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
    nlParams = PyTrilinos.NOX.Epetra.defaultNonlinearParameters(comm,
                                                                verbosity,
                                                                outputPrec,
                                                                maxIterations,
                                                                tolerance)
    direction    = nlParams["Direction"]
    newton       = nlParams["Newton"]
    linearSolver = nlParams["Linear Solver"]
    linearSolver["Output Frequency"           ] = 1
    linearSolver["Preconditioner"             ] = "None"
    linearSolver["Preconditioner Operator"    ] = "Use Jacobian"
    linearSolver["Size of Krylov Subspace"    ] = 100
    linearSolver["Tolerance"                  ] = 1e-08
    linearSolver["Zero Initial Guess"         ] = False
    linearSolver["Compute Scaling Manually"   ] = True
    linearSolver["Throw Error on Prec Failure"] = True
    linearSolver["RCM Reordering"             ] = "Disabled"
    linearSolver["Orthogonalization"          ] = "Classical"
    linearSolver["Convergence Test"           ] = "r0"
    linearSolver["Preconditioner Reuse Policy"] = "Rebuild"
    newton["Linear Solver"] = linearSolver
    direction["Newton"] = newton
    noxParams = {"Tolerance"        : tolerance,
                 "Printing"         : nlParams["Printing"],
                 "Nonlinear Solver" : nlParams["Nonlinear Solver"],
                 "Direction"        : direction
                 }

    predictor   = {"Method" : "Secant"}
    bifurcation = {"Type"   : "None"  }
    stepSize    = {"Method"            : "Adaptive",  
                   "Initial Step Size" : 0.01,
                   "Min Step Size"     : 1.0e-3,
                   "Max Step Size"     : .02
                   }
    eigensolver = {"Method"                : "Anasazi",
                   "Sorting Order"         : "LM",
                   "Block Size"            : 4,
                   "Num Blocks"            : 100,
                   "Num Eigenvalues"       : 15,
                   "Convergence Tolerance" : 1e-11,
                   "Step Size"             : 20,
                   "Maximum Restarts"      : 20,
                   "Maximum Iterations"    : 500,
                   "Operator"              : "Jacobian Inverse",
                   "Operator"              : "Shift-Invert",
                   "Shift"                 : 0.03,
                   "Cayley Pole"           : 0.001,
                   "Cayley Zero"           : 0.02,
                   "Symmetric"             : False
                   }
    stepper = {"Continuation Method"      : "Arc Length",
               "Continuation Parameter"   : "sigma",
               "Initial Value"            : 0.0,
               "Max Value"                : 10.0,
               "Min Value"                : 0.0,
               "Max Nonlinear Iterations" : 100,
               "Max Steps"                : 100,
               "Compute Eigenvalues"      : False,
               "Eigensolver"              : eigensolver
               }

    locaParams = {"Enable Arc Length Scaling"               : False,
                  "Goal Arc Length Parameter Contribution"  : 0.8,
                  "Max Arc Length Parameter Contribution"   : 0.8,
                  "Initial Scale Factor"                    : 1.0,
                  "Min Scale Factor"                        : 1e-7,
                  "Enable Tangent Factor Step Size Scaling" : True,
                  "Min Tangent Factor"                      : 0.8,
                  "Tangent Factor Exponent"                 : 1.5,
                  "Predictor"                               : predictor,
                  "Bifurcation"                             : bifurcation, 
                  "Step Size"                               : stepSize,
                  "Stepper"                                 : stepper
                  }

    result = {"NOX"  : noxParams,
              "LOCA" : locaParams
              }

    return result
%}
