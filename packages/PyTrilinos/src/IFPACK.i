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
PyTrilinos.IFPACK is the python interface to the Trilinos
preconditioner package IFPACK:

    http://trilinos.sandia.gov/packages/ifpack

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
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"

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
#endif

// Teuchos Python utility code
#ifdef HAVE_TEUCHOS
#include "Teuchos_PythonParameter.h"
#endif

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

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include IFPACK documentation
%include "IFPACK_dox.i"

// External Trilinos modules
#ifdef HAVE_TEUCHOS
%import "Teuchos.i"
#endif
#ifdef HAVE_EPETRA
%ignore Epetra_Version();
%import "Epetra.i"
#endif

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unkown C++ exception");
  }
}

////////////////////////
// I/O Stream support //
////////////////////////
// #pragma SWIG nowarn=401
//%include "std_iostream.i"

//////////////////////////////////
// IFPACK configuration support //
//////////////////////////////////
#undef  PACKAGE_BUGREPORT
%ignore PACKAGE_BUGREPORT;
#undef  PACKAGE_NAME
%ignore PACKAGE_NAME;
#undef  PACKAGE_STRING
%ignore PACKAGE_STRING;
#undef  PACKAGE_TARNAME
%ignore PACKAGE_TARNAME;
#undef  PACKAGE_VERSION
%ignore PACKAGE_VERSION;
%include "Ifpack_config.h"
%include "Ifpack_ConfigDefs.h"

////////////////////////////
// IFPACK_Version support //
////////////////////////////
%include "Ifpack_Version.h"
%pythoncode %{
Version = Ifpack_Version
__version__ = Version().split()[3]
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
%teuchos_rcp(Ifpack_Preconditioner)
%ignore operator<<(ostream &, const Ifpack_Preconditioner &);
%ignore Ifpack_Preconditioner::Condest() const;
%rename(Preconditioner) Ifpack_Preconditioner;
%include "Ifpack_Preconditioner.h"
%extend Ifpack_Preconditioner
{
  std::string __str__()
  {
    std::stringstream os;
    os << *self;
    return os.str();
  }
  //void __del__()
  //{
  //  delete self;
  //}
}

///////////////////////
// IFPACK_IC support //
///////////////////////
%teuchos_rcp(Ifpack_IC)
%ignore Ifpack_IC::Condest() const;
%rename(IC) Ifpack_IC;
%include "Ifpack_IC.h"

////////////////////////
// IFPACK_ICT support //
////////////////////////
%teuchos_rcp(Ifpack_ICT)
%ignore Ifpack_ICT::Condest() const;
%rename(ICT) Ifpack_ICT;
%include "Ifpack_ICT.h"

////////////////////////
// IFPACK_ILU support //
////////////////////////
%teuchos_rcp(Ifpack_ILU)
%ignore Ifpack_ILU::Condest() const;
%rename(ILU) Ifpack_ILU;
%include "Ifpack_ILU.h"

/////////////////////////
// IFPACK_ILUT support //
/////////////////////////
%teuchos_rcp(Ifpack_ILUT)
%ignore Ifpack_ILUT::Condest() const;
%rename(ILUT) Ifpack_ILUT;
%include "Ifpack_ILUT.h"

////////////////////////////////////
// IFPACK_PointRelaxation support //
////////////////////////////////////
%teuchos_rcp(Ifpack_PointRelaxation)
%ignore Ifpack_PointRelaxation::Condest() const;
%rename(PointRelaxation) Ifpack_PointRelaxation;
%include "Ifpack_PointRelaxation.h"

///////////////////////////
// IFPACK_Amesos support //
///////////////////////////
%teuchos_rcp(Ifpack_Amesos)
%ignore Ifpack_Amesos::Condest() const;
%rename(Amesos) Ifpack_Amesos;
%include "Ifpack_Amesos.h"

////////////////////////////
// IFPACK factory support //
////////////////////////////
%ignore operator>>(std::istream &, Ifpack::EPrecType &);
%newobject Ifpack::Create;
%rename(Factory) Ifpack;
%include "Ifpack.h"

// Turn off the exception handling
%exception;
