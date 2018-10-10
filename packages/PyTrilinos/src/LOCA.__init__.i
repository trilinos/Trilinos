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

%define %loca_docstring
"
PyTrilinos.LOCA is the python interface to the Trilinos continuation
algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA is to provide a library of continuation
algorithms.  It include files the following sub-modules:

    * Abstract           - Abstract continuation problem base classes
    * Extended           - Classes that extend NOX.Abstract classes to
                           handle an arbitrary number of multi-vectors
                           and scalars
    * MultiContinuation  - Groups and vectors for multi-parameter continuation
    * TimeDependent      - Abstract group for time dependent problems with a
                           mass matrix
    * TurningPoint       - Groups and vectors for turning point bifurcations
    * Hopf               - Groups and vectors for Hopf bifurcations
    * Pitchfork          - Groups and vectors for pitchfork bifurcations
    * Homotopy           - Groups that allow for Homotopy to be applied
    * PhaseTransition    - Groups and vectors for phase transition bifurcations
    * Parameter          - Centralized library for setting/retrieving numerical
                           parameter values in application codes
    * BorderedSolver     - Strategies for solving bordered systems of equations
    * BorderedSystem     - Interface for groups that are bordered systems
    * Bifurcation        - Strategies for creating bifurcation objects
    * StatusTest         - Status checkers
    * StepSize           - Collection of step size control strategies
    * MultiPredictor     - Predictor direction strategies

and classes:

    * GlobalData      - Container class that holds ref-count pointers to
                        'global' objects, i.e., objects that nearly every LOCA
                        object will need access to
    * ErrorCheck      - Error checking algorithm for NOX/LOCA routines
    * Factory         - Provides a single location for instantiating various
                        strategies based on parameter list choices
    * DerivUtils      - Generic derivative computation class to compute various
                        derivatives via finite differencing
    * Stepper         - Implementation of LOCA.Abstract.Iterator for computing
                        points along a continuation curve
    * ParameterVector - LOCA's container for holding a set of parameters that
                        are used by the LOCA continuation routines
"
%enddef

%module(package      = "PyTrilinos.LOCA",
        directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_docstring) __init__

%{
// System include files
#include <sstream>

// PyTrilinos include files
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_LinearProblem.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"
#endif

// LOCA include files
#include "PyTrilinos_LOCA_Headers.hpp"
#include "PyTrilinos_LOCA_Hopf_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

%pythoncode
%{
# Fix ___init__ ambiguity
del ___init__
from . import ___init__
%}

// Ignore/renames
%ignore *::operator=;
%ignore *::operator[];
%ignore operator<<(ostream&, const LOCA::ParameterVector&);
%rename(Print) LOCA::ParameterVector::print(ostream& stream) const;

// SWIG library include files
%include "stl.i"

// Trilinos interface import.  Note: Teuchos.i turns off warnings for
// nested classes, so we do not have to do it again.
%import "Teuchos.i"

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

// NOX interface file imports.
%import "NOX.Abstract.i"
%import "NOX.StatusTest.i"

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
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Import NOX and LOCA sub-modules
%pythoncode
%{
import PyTrilinos.NOX

# LOCA sub-modules
__all__ = ['Extended',
           'MultiContinuation',
           'TimeDependent',
           'TurningPoint',
           'Hopf',
           'Pitchfork',
           'Homotopy',
           'PhaseTransition',
           'Abstract',
           'Parameter',
           'BorderedSolver',
           'BorderedSystem',
           'Bifurcation',
           'StatusTest',
           'StepSize',
           'MultiPredictor',
           'Eigensolver',
           'EigenvalueSort',
           'SaveEigenData',
           'AnasaziOperator'
           ]
from . import Extended
from . import MultiContinuation
from . import TimeDependent
from . import TurningPoint
from . import Hopf
from . import Pitchfork
from . import Homotopy
from . import PhaseTransition
from . import Parameter
from . import BorderedSolver
from . import BorderedSystem
from . import Bifurcation
from . import StatusTest
from . import StepSize
from . import MultiPredictor
from . import Eigensolver
from . import EigenvalueSort
from . import SaveEigenData
from . import AnasaziOperator
%}

// Techos::RCP handling
%teuchos_rcp(LOCA::GlobalData)
%teuchos_rcp(LOCA::ErrorCheck)
%teuchos_rcp(LOCA::Factory)
%teuchos_rcp(LOCA::Stepper)
%teuchos_rcp(LOCA::DerivUtils)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)

// LOCA GlobalData class
%include "LOCA_GlobalData.H"

// LOCA ErrorCheck class
%include "LOCA_ErrorCheck.H"

// LOCA Factory class
%include "LOCA_Factory.H"

// LOCA DerivUtils class
%include "LOCA_DerivUtils.H"

// The LOCA::Stepper class derives from LOCA::Abstract::Iterator, so
// import it here
%teuchos_rcp(LOCA::Abstract::Iterator)
%import(module="Abstract") "LOCA_Abstract_Iterator.H"

// At this point, 'Abstract' might be 'NOX.Abstract' (depending on the
// Python version and related import rules), but we need it to be
// 'LOCA.Abstract'
%pythoncode
%{
import os.path
if 'NOX' in Abstract.__file__.split(os.path.sep):
  del Abstract
  from . import Abstract
%}

// LOCA Stepper class
%teuchos_rcp(LOCA::Stepper)
%feature("director") LOCA::Stepper;
// Ignore the deprecated constructor
// %ignore LOCA::Stepper::Stepper(const Teuchos::RCP< LOCA::GlobalData >&,
//                                const Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >&,
//                                const Teuchos::RCP< NOX::StatusTest::Generic >&,
//                                const Teuchos::RCP< Teuchos::ParameterList >&);
%include "LOCA_Stepper.H"

// LOCA ParameterVector class
%include "LOCA_Parameter_Vector.H"

// LOCA.Epetra
#ifdef HAVE_PYTRILINOS_NOX_EPETRA
%pythoncode
%{

# Epetra namespace
__all__.append("Epetra")
from . import Epetra
%}
#endif
