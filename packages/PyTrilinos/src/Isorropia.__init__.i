// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2010) Sandia Corporation
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

%define %isorropia_docstring
"
PyTrilinos.Isorropia is the python interface to the Trilinos
partitioning and load balancing package Isorropia:

    http://trilinos.sandia.gov/packages/isorropia

The purpose of Isorropia is to ....
"
%enddef

%module(package   = "PyTrilinos.Isorropia",
	autodoc   = "1",
	docstring = %isorropia_docstring) __init__

%{
// PyTrilinos configuration
#include "PyTrilinos_config.h"

// Isorropia includes
#include "Isorropia_Version.hpp"
#include "Isorropia_Operator.hpp"
#include "Isorropia_Colorer.hpp"
#include "Isorropia_Partitioner.hpp"
#include "Isorropia_Partitioner2D.hpp"
#include "Isorropia_Redistributor.hpp"
#include "Isorropia_CostDescriber.hpp"
#include "Isorropia_Orderer.hpp"
#include "Isorropia_LevelScheduler.hpp"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "PyTrilinos_Teuchos_Util.h"
%}

// General ignore directives
%ignore operator<<;
%ignore *::operator=;
%ignore *::operator[];

// Auto-documentation feature
%feature("autodoc", "1");

// Include Isorropia documentation (this file will need to be
// generated before it can be included)
%include "Isorropia_dox.i"

// Trilinos interface import
%import "Teuchos.i"

///////////////////////////////
// Isorropia Version support //
///////////////////////////////
%include "Isorropia_Version.hpp"
%pythoncode
%{
__version__ = Isorropia_Version().split()[3]
%}

/////////////////////////////////
// Isorropia::Operator support //
/////////////////////////////////
%teuchos_rcp(Isorropia::Operator)
%include "Isorropia_Operator.hpp"

////////////////////////////////
// Isorropia::Colorer support //
////////////////////////////////
%teuchos_rcp(Isorropia::Colorer)
%extend Isorropia::Colorer
{
  PyObject * elemsWithColor(int color)
  {
    int length = self->numElemsWithColor(color);
    npy_intp dims[1] = { length };
    PyObject * elemArray = PyArray_SimpleNew(1, dims, NPY_INT);
    if (PyErr_Occurred()) return NULL;
    int * elementList = (int*) array_data(elemArray);
    self->elemsWithColor(color, elementList, length);
    return elemArray;
  }
}
%ignore Isorropia::Colorer::elemsWithColor;
%include "Isorropia_Colorer.hpp"

////////////////////////////////////
// Isorropia::Partitioner support //
////////////////////////////////////
%teuchos_rcp(Isorropia::Partitioner)
%include "Isorropia_Partitioner.hpp"

//////////////////////////////////////
// Isorropia::Partitioner2D support //
//////////////////////////////////////
//%teuchos_rcp(Isorropia::Partitioner2D)
//%include "Isorropia_Partitioner2D.hpp"

//////////////////////////////////////
// Isorropia::Redistributor support //
//////////////////////////////////////
//%teuchos_rcp(Isorropia::Redistributor)
%include "Isorropia_Redistributor.hpp"

//////////////////////////////////////
// Isorropia::CostDescriber support //
//////////////////////////////////////
%teuchos_rcp(Isorropia::CostDescriber)
%include "Isorropia_CostDescriber.hpp"

////////////////////////////////
// Isorropia::Orderer support //
////////////////////////////////
%teuchos_rcp(Isorropia::Orderer)
%include "Isorropia_Orderer.hpp"

///////////////////////////////////////
// Isorropia::LevelScheduler support //
///////////////////////////////////////
%teuchos_rcp(Isorropia::LevelScheduler)
%include "Isorropia_LevelScheduler.hpp"

// Isorropia namespace imports
%pythoncode
%{
# Epetra namespace
__all__ = ['Epetra']
import NestedEpetra as Epetra
%}
