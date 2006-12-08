// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1") __init__

%{
// Teuchos includes
#include "Teuchos_PythonParameter.hpp"

// Epetra includes
#include "Epetra_BLAS.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_JadMatrix.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

// Local includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"

// Namespace flattening
using namespace NOX;
using namespace NOX::Abstract;
using namespace NOX::Epetra;
%}

// The following #define is to change the name of a NOX method
// argument that conflicts with a SWIG director method argument
#define result nox_result

// Ignore directives
%ignore *::print(std::ostream &) const;
%ignore *::print(std::ostream &, int) const;
%ignore *::operator=;
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);
%ignore NOX::Epetra::Vector::getEpetraVector() const;

// Rename directives
%rename(Group_None) NOX::Epetra::Group::None;

// Feature directives
%feature("director") NOX::Epetra::LinearSystem;

// SWIG library includes
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Support for Teuchos::RefCountPtrs
TEUCHOS_RCP_TYPEMAPS(NOX::Epetra::LinearSystem)
TEUCHOS_RCP_TYPEMAPS(Epetra_Operator)

// Make Epetra_Vector and NOX::Epetra::Vector input arguments
// interchangeable
%typemap(in) NOX::Epetra::Vector & (void* argp=0, int res=0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor(Epetra_Vector*), %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp) { %argument_nullref("$type", $symname, $argnum); }
    $1 = new NOX::Epetra::Vector(*%reinterpret_cast( argp, Epetra_Vector*));
  } else {
    $1 = %reinterpret_cast(argp, NOX::Epetra::Vector*);
  }
}

%typecheck(1190) NOX::Epetra::Vector & {
  static void * argp = 0;
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor(Epetra_Vector*),
                                         %convertptr_flags)) ? 1 : 0;
}

%typemap(freearg) NOX::Epetra::Vector {
  delete $1;
}

// Epetra imports
%import  "Epetra_SrcDistObject.h"
%import  "Epetra_Operator.h"
%import  "Epetra_RowMatrix.h"

// NOX imports
%import "NOX_Abstract_Group.H"
%import "NOX_Abstract_Vector.H"

// NOX::Epetra::Interface imports
%import "NOX.Epetra.Interface.i"

// NOX::Epetra includes
using namespace std;
%include "NOX_Epetra_Group.H"
%include "NOX_Epetra_Vector.H"
%include "NOX_Epetra_FiniteDifference.H"
%include "NOX_Epetra_FiniteDifferenceColoring.H"
%include "NOX_Epetra_MatrixFree.H"
%include "NOX_Epetra_LinearSystem.H"
%include "NOX_Epetra_LinearSystem_AztecOO.H"
