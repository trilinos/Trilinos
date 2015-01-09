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

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

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
#include "numpy_include.hpp"
#include "PyTrilinos_Teuchos_Util.hpp"
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
%teuchos_rcp(Isorropia::Redistributor)
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
import IsorropiaEpetra as Epetra
%}
