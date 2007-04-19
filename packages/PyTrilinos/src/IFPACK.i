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

%define %ifpack_docstring
"
PyTrilinos.IFPACK is the python interface to Trilinos package IFPACK:

    http://software.sandia.gov/trilinos/packages/ifpack

The purpose of IFPACK is to provide incomplete foctorization
preconditioners to Trilinos.  Note that the C++ version of IFPACK uses
the prefix 'Ifpack_' which has been stripped from the python version.

IFPACK provides the following user-level classes:

    * Factory         - A factory for creating IFPACK preconditioners
    * Preconditioner  - Pure virtual base class for defining interface
    * IC              - Incomplete Cholesky preconditioner
    * ICT             - Incomplete Cholesky preconditioner w/threshold
    * ILU             - Incomplete lower/upper preconditioner
    * ILUT            - Incomplete lower/upper preconditioner w/threshold
    * PointRelaxation - Point relaxation predonditioner
    * Amesos          - Use Amesos factorizations as preconditioners

and functions:

    * AnalyzeMatrix          - Analyze the basic properties of a matrix
    * AnalyzeMatrixElements  - Analyze the distribution of values of a matrix
    * AnalyzeVectorElements  - Analyze the distribution of values of a vector
    * PrintSparsity          - Create PS file with sparsity pattern of matrix

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exIFPACK.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %ifpack_docstring) IFPACK

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_JadOperator.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"

// Epetra python includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

// Teuchos Python utility code
#include "Teuchos_PythonParameter.h"

// IFPACK includes
#include "Ifpack.h"
#include "Ifpack_Version.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Amesos.h"

%}

// Auto-documentation feature
%feature("autodoc", "1");

// External Trilinos modules
using namespace std;
%import "Teuchos.i"
%ignore Epetra_Version();
%import "Epetra.i"

//////////////////////////////////
// IFPACK configuration support //
//////////////////////////////////
%include "Ifpack_config.h"
%include "Ifpack_ConfigDefs.h"

////////////////////////////
// IFPACK factory support //
////////////////////////////
%rename(Factory) Ifpack;
%include "Ifpack.h"

////////////////////////////
// IFPACK_Version support //
////////////////////////////
%rename(Version) Ifpack_Version;
%include "Ifpack_Version.h"
%pythoncode %{
__version__ = Version().split()[2]
%}

//////////////////////////
// IFPACK_Utils support //
//////////////////////////
%rename(AnalyzeMatrix        ) Ifpack_Analyze;
%rename(AnalyzeMatrixElements) Ifpack_AnalyzeMatrixElements;
%rename(AnalyzeVectorElements) Ifpack_AnalyzeVectorElements;
%rename(PrintSparsity        ) Ifpack_PrintSparsity;
%include "Ifpack_Utils.h"

///////////////////////////////////
// IFPACK_Preconditioner support //
///////////////////////////////////
%ignore operator<<(ostream &, const Ifpack_Preconditioner &);
%ignore Ifpack_Preconditioner::Condest() const;
%rename(Preconditioner) Ifpack_Preconditioner;
%include "Ifpack_Preconditioner.h"
%extend Ifpack_Preconditioner {
  string __str__() {
    stringstream os;
    os << *self;
    return os.str();
  }
  void __del__()
  {
    delete self;
  }
}

///////////////////////
// IFPACK_IC support //
///////////////////////
%ignore Ifpack_IC::Condest() const;
%rename(IC) Ifpack_IC;
%include "Ifpack_IC.h"

////////////////////////
// IFPACK_ICT support //
////////////////////////
%ignore Ifpack_ICT::Condest() const;
%rename(ICT) Ifpack_ICT;
%include "Ifpack_ICT.h"

////////////////////////
// IFPACK_ILU support //
////////////////////////
%ignore Ifpack_ILU::Condest() const;
%rename(ILU) Ifpack_ILU;
%include "Ifpack_ILU.h"

/////////////////////////
// IFPACK_ILUT support //
/////////////////////////
%ignore Ifpack_ILUT::Condest() const;
%rename(ILUT) Ifpack_ILUT;
%include "Ifpack_ILUT.h"

////////////////////////////////////
// IFPACK_PointRelaxation support //
////////////////////////////////////
%ignore Ifpack_PointRelaxation::Condest() const;
%rename(PointRelaxation) Ifpack_PointRelaxation;
%include "Ifpack_PointRelaxation.h"

///////////////////////////
// IFPACK_Amesos support //
///////////////////////////
%ignore Ifpack_Amesos::Condest() const;
%rename(Amesos) Ifpack_Amesos;
%include "Ifpack_Amesos.h"
