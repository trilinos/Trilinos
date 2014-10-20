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

%define %pliris_docstring
"
PyTrilinos.Pliris is the python interface to the Trilinos package
Pliris, an LU solver for dense matrices:

    http://trilinos.sandia.gov/packages/pliris

The purpose of Pliris is to provide an object-oriented interface to an
LU solver for dense matrices on parallel platforms. These matrices are
double precision real matrices distributed on a parallel machine.  The
python version of the Pliris package supports the following class:

    * Pliris  - Primary solver class.
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %pliris_docstring) Pliris

%{
// System includes
// #include <iostream>
// #include <sstream>
// #include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#ifdef HAVE_EPETRA
#undef HAVE_INTTYPES_H
#undef HAVE_STDINT_H
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
// #include "Epetra_MultiVector.h"
// #include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"

// Epetra python includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"

// Pliris includes
// #include "Pliris_config.h"
#include "Pliris.h"

#endif

%}

// Auto-documentation feature
%feature("autodoc", "1");

// Include Pliris documentation
// %include "Pliris_dox.i"

// External Trilinos modules
#ifdef HAVE_EPETRA
%ignore Epetra_Version();
%import "Epetra.i"
#endif

///////////////////////////////////
// Pliris configuration support //
///////////////////////////////////
%include "Pliris_config.h"
%pythoncode
{
__version__ = PACKAGE_VERSION
}

////////////////////
// Pliris support //
////////////////////
%include "Pliris.h"
