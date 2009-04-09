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

// Teuchos include
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_Version.H"
#include "NOX_Utils.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_LineSearchBased.H"
#include "NOX_Solver_TrustRegionBased.H"
#include "NOX_Solver_InexactTrustRegionBased.H"
#include "NOX_Solver_TensorBased.H"
#include "NOX_Solver_Factory.H"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"

// Local includes
#include "NumPyImporter.h"
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

// Trilinos interface file imports.
// The Teuchos.py file that this %import will try to import in python
// will be one directory up from this python module.  So we need to
// add the parent directory to the search path.
%pythoncode
{
import os.path, sys
currentDir,dummy = os.path.split(__file__)
sys.path.append(os.path.normpath(os.path.join(currentDir,"..")))
}
%import "Teuchos.i"
// Note: Teuchos.i turns off warnings for nested classes, so we do not
// have to do it again.

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
%import "NOX.Abstract.i"
%import "NOX.Solver.i"
%import "NOX.StatusTest.i"

// Python code for the NOX module
%pythoncode
%{
try:
    import Epetra
except ImportError, e:
    if str(e) != "No module named Epetra":
        print e
%}
