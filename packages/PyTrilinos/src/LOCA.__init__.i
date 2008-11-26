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

%define %loca_docstring
"
PyTrilinos.LOCA is the python interface to the Trilinos continuation
algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA is to provide a library of continuation
algorithms.  This module is not currently supported, but the plan is
to reactivate it soon.
"
%enddef

%module(package   = "PyTrilinos.LOCA",
	docstring = %loca_docstring) __init__

%{
// System includes
#include <sstream>

// Teuchos include
#include "Teuchos_PythonParameter.h"

// LOCA includes
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"

// Local includes
#include "NumPyImporter.h"
%}

// Ignore/renames
%ignore *::operator=;
%ignore *::operator[];
%ignore operator<<(ostream&, const LOCA::ParameterVector&);
%rename(Print) LOCA::ParameterVector::print(ostream& stream) const;

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

// LOCA interface includes
%include "LOCA_Abstract_Iterator.H"
%include "LOCA_Stepper.H"
%include "LOCA_Parameter_Vector.H"
