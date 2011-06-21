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

%define %komplex_docstring
"
PyTrilinos.Komplex is the python interface to the Trilinos complex
linear algebra package Komplex:

    http://trilinos.sandia.gov/packages/komplex

The purpose of Komplex is to define complex (real + imaginary) linear
algebra problems using real-valued Epetra vectors and matrix
operators.  The python version of the Komplex package supports the
following class:

    * LinearProblem - Encapsulate all elements of a complex linear
                      algebra problem
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %komplex_docstring) Komplex

%{
// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_LinearProblem.h"
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

// Epetra python includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"

// Komplex includes
#include "Komplex_Version.h"
#include "Komplex_LinearProblem.h"

#endif

%}

// Auto-documentation feature
%feature("autodoc", "1");

// Include Komplex documentation
%include "Komplex_dox.i"

// SWIG library includes
using std::string;
%include "stl.i"

///////////////////////////////////
// Komplex configuration support //
///////////////////////////////////
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
%include "Komplex_config.h"
%include "PyTrilinos_config.h"

// Teuchos::RCP<> support
#ifdef TEUCHOS
%include "Teuchos_RCP.i"
#endif

// External Trilinos modules
#ifdef HAVE_EPETRA
%ignore Epetra_Version();
%import "Epetra.i"
#endif

/////////////////////////////
// Komplex Version support //
/////////////////////////////
%rename(Version) Komplex_Version;
%include "Komplex_Version.h"
%pythoncode
{
__version__ = Version().split()[2]
}

///////////////////////////////////
// Komplex LinearProblem support //
///////////////////////////////////
%include "Komplex_LinearProblem.h"
// I don't use %rename here, because it conflicts with
// Epetra_LinearProblem.  I just create a new symbol within the
// Komplex namespace
%pythoncode
{
LinearProblem = Komplex_LinearProblem
}
