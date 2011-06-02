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

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_epetra_docstring) __init__

%{
// System includes
#include <vector>

// Configuration
#include "PyTrilinos_config.h"

// Teuchos includes
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_BLAS.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_DLLExportMacro.h"

// EpetraExt includes
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_ModelEvaluator.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_Scaling.H"
#include "NOX_Epetra_LinearSystem.H"
#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_ModelEvaluatorInterface.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyVector.h"

// Namespace flattening
using Teuchos::RCP;
using Teuchos::rcp;
using namespace NOX;
using namespace NOX::Abstract;
using namespace NOX::Epetra;
%}

// Configuration
%include "Epetra_DLLExportMacro.h"

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

// SWIG library includes
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"

// Support for Teuchos::RCPs
%teuchos_rcp(NOX::Epetra::LinearSystem)
%teuchos_rcp(NOX::Epetra::LinearSystemAztecOO)
%teuchos_rcp(NOX::Epetra::Scaling)
%teuchos_rcp(NOX::Epetra::VectorSpace)

//////////////
// Typemaps //
//////////////

// Make Epetra_Vector and NOX::Epetra::Vector input arguments
// interchangeable
%typemap(in) NOX::Epetra::Vector &
(void* argp=0, int res=0, Teuchos::RCP< Epetra_NumPyVector > tempshared, bool cleanup=false)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res))
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(Teuchos::RCP< Epetra_NumPyVector > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp)
    {
      %argument_nullref("$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< Epetra_NumPyVector > *);
      delete %reinterpret_cast(argp, Teuchos::RCP< Epetra_NumPyVector > *);
      $1 = new NOX::Epetra::Vector(Teuchos::rcp_dynamic_cast< Epetra_Vector >(tempshared),
				   NOX::Epetra::Vector::CreateView);
      cleanup = true;
    }
    else
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< Epetra_NumPyVector > *);
      $1 = new NOX::Epetra::Vector(Teuchos::rcp_dynamic_cast< Epetra_Vector >(tempshared),
				   NOX::Epetra::Vector::CreateView);
      cleanup = true;
    }
  }
  else
  {
    $1 = %reinterpret_cast(argp, NOX::Epetra::Vector*);
  }
}
%typecheck(1190) NOX::Epetra::Vector &
{
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, 0, $descriptor, 0)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtrAndOwn($input, 0,
					       $descriptor(Teuchos::RCP< Epetra_NumPyVector > *),
					       %convertptr_flags, 0)) ? 1 : 0;
}
%typemap(freearg) NOX::Epetra::Vector &
{
  if (cleanup$argnum) delete $1;
}

// Convert NOX::Abstract::Vector return arguments to Epetra.Vectors
%typemap(out) NOX::Abstract::Vector &
(NOX::Epetra::Vector* nevResult  = NULL)
{
  nevResult = dynamic_cast<NOX::Epetra::Vector*>($1);
  if (nevResult == NULL)
  {
    // If we cannot downcast, then return the NOX::Abstract::Vector
    $result = SWIG_NewPointerObj((void*)&$1, $descriptor, 1);
  }
  else
  {
    Epetra_NumPyVector enpvResult(View, nevResult->getEpetraVector(), 0);
    Teuchos::RCP< Epetra_NumPyVector > *smartresult = 
      new Teuchos::RCP< Epetra_NumPyVector >(enpvResult, bool($owner));
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< Epetra_NumPyVector > *),
				   SWIG_POINTER_OWN));
  }
}

// Convert Epetra_Vector return arguments to Epetra.Vectors (now provided by Epetra_Base.i)
// %typemap(out) Epetra_Vector & (Epetra_NumPyVector* enpvResult = NULL)
// {
//   enpvResult = new Epetra_NumPyVector(View, *$1, 0);
//   $result = SWIG_NewPointerObj((void*)enpvResult, $descriptor(Epetra_NumPyVector*), 1);
// }

// Convert NOX::Epetra::LinearSystem objects to
// NOX::Epetra::LinearSystemAztecOO
%typemap(out) Teuchos::RCP< NOX::Epetra::LinearSystem >
(NOX::Epetra::LinearSystem*        nelsPtr     = NULL,
 NOX::Epetra::LinearSystemAztecOO* nelsaResult = NULL)
{
  nelsPtr = $1.get();
  nelsaResult = dynamic_cast< NOX::Epetra::LinearSystemAztecOO*>(nelsPtr);
  if (nelsaResult == NULL)
  {
    //If we cannot downcast then return the NOX::Epetra::LinearSystem
    %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystem > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *smartresult =
      new Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO >(nelsaResult);
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *),
				   SWIG_POINTER_OWN));
  }
}

%typemap(out) Teuchos::RCP< const NOX::Epetra::LinearSystem >
(const NOX::Epetra::LinearSystem*        nelsPtr     = NULL,
 const NOX::Epetra::LinearSystemAztecOO* nelsaResult = NULL)
{
  nelsPtr = $1.get();
  nelsaResult = dynamic_cast< const NOX::Epetra::LinearSystemAztecOO*>(nelsPtr);
  if (nelsaResult == NULL)
  {
    //If we cannot downcast then return the NOX::Epetra::LinearSystem
    %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystem > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    Teuchos::RCP< const NOX::Epetra::LinearSystemAztecOO > *smartresult =
      new Teuchos::RCP< const NOX::Epetra::LinearSystemAztecOO >(nelsaResult);
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *),
				   SWIG_POINTER_OWN));
  }
}

// Epetra includes.  This is a potential source of problems.  The
// simple thing to do is to add an "%import 'Epetra.i'" here.  If I do
// that, strange things start to happen: other, seemingly unrelated
// wrappers start to seg fault.  I do not have a good explanation for
// it.  By %include-ing the following, I bypass inserting an "import
// PyTrilinos.Epetra" into the resulting .py file (because Epetra_*.i
// files do not have a %module directive).  This seems to provide the
// functionality I need without causing whatever confusion is at risk
// here.
%include "Epetra_Base.i"
%teuchos_rcp_epetra_numpy(IntSerialDenseVector)
%teuchos_rcp_epetra_numpy(SerialDenseVector)
%teuchos_rcp_epetra_numpy(SerialDenseMatrix)
%include "Epetra_Comm.i"
%include "Epetra_Maps.i"
%include "Epetra_Dist.i"
%teuchos_rcp_epetra_numpy(MultiVector)
%teuchos_rcp_epetra_numpy(Vector)
%include "Epetra_Graphs.i"
%include "Epetra_Operators.i"

// EpetraExt import
%ignore EpetraExt::Add;
%include "EpetraExt.i"

// NOX import
%import "NOX.Abstract.i"

// NOX::Epetra::Interface imports
%import "NOX.Epetra.Interface.i"

//////////////////////////////
// NOX.Epetra.Group support //
//////////////////////////////
%teuchos_rcp(NOX::Epetra::Group)
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
%teuchos_rcp(NOX::Epetra::FiniteDifference)
%include "NOX_Epetra_FiniteDifference.H"

/////////////////////////////////////////////////
// NOX.Epetra.FiniteDifferenceColoring support //
/////////////////////////////////////////////////
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

///////////////////////////////////
// NOX.Epetra.MatrixFree support //
///////////////////////////////////
%teuchos_rcp(NOX::Epetra::MatrixFree)
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

// Turn off the exception handling
%exception;

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
        assert isinstance(jacobian    , (PyTrilinos.Epetra.Operator, Operator))
    if precInterface is not None:
        assert isinstance(precInterface , Interface.Preconditioner  )
        assert isinstance(preconditioner, (PyTrilinos.Epetra.Operator, Operator))

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
        assert isinstance(jacobian    , (PyTrilinos.Epetra.Operator, Operator))
    if precInterface is not None:
        assert isinstance(precInterface , Interface.Preconditioner  )
        assert isinstance(preconditioner, (PyTrilinos.Epetra.Operator, Operator))

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
