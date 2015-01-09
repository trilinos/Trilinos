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
#include "numpy_include.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"

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
