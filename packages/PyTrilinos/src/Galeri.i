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

// This documentation string will be the python help facility help
// string
%define %galeri_docstring
"Galeri: Matrix Generation Package.
The Galeri module allows an easy creation of Epetra_Map's and
Epetra_CrsMatrix's.  Use the python help() facility for local documentation on
classes and methods, or see the on-line documentation for more in-depth
information. Also give a look to the examples in galeri/python/example.

The most important classes of the Galeri module are:
- Galeri.CreateMap()
- Galeri.CreateCrsMatrix()
- Galeri.ReadHB()

Galeri requires the Epetra and Teuchos modules of PyTrilinos.
"

%enddef

// Define the module name, its package and documentation string
%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %galeri_docstring) Galeri

%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"

// Epetra Python utility code
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Galeri includes
#include "Galeri_Version.h"
#include "Galeri_Utils.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_ReadHB.h"
%}

// Trilinos package imports
%import "Teuchos.i"
%import "Epetra.i"

// Turn on autodocumentation
%feature("autodoc", "1");

// Typemap support for STL
using namespace std;
%include "stl.i"

////////////////////////////
// Galeri_Version support //
////////////////////////////
%include "Galeri_Version.h"
%pythoncode %{
__version__ = Galeri_Version().split()[2]
%}

//////////////////////////
// Galeri_Utils support //
//////////////////////////
%include "Galeri_Utils.h"

/////////////////////////
// Galeri_Maps support //
/////////////////////////
%include "Galeri_Maps.h"

////////////////////////////////
// Galeri_CrsMatrices support //
////////////////////////////////
%include "Galeri_CrsMatrices.h"

////////////////////////////////
// Galeri_VbrMatrices support //
////////////////////////////////
%include "Galeri_VbrMatrices.h"

///////////////////////////
// Galeri_ReadHB support //
///////////////////////////
%include "Galeri_ReadHB.h"
