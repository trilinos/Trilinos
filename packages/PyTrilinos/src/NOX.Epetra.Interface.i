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

%define %nox_epetra_interface_docstring
"
PyTrilinos.NOX.Epetra.Interface is the python interface to the
Epetra::Interface namespace of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Epetra.Interface is to provide base classes the
user should derive from in order to define the nonlinear function to
be solved, and if needed, its Jacobian and the desired preconditioner.

NOX.Epetra.Interface provides the following user-level classes:

    * Required        - Required class for computing the nonlinear function
    * Jacobian        - Class for computing the Jacobian (if needed)
    * Preconditioner  - Class for computing the preconditioner (if needed)
"
%enddef

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_epetra_interface_docstring) Interface

%{
// NumPy includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
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

// Local Epetra includes
#include "Epetra_PyUtil.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"

// NOX::Epetra::Interface includes
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
%}

// General ignore directives
// %ignore *::operator=;   // temp removal

// Include NOX documentation
%include "NOX_dox.i"

// STL support
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Epetra module imports
%import "Epetra.i"

// Exception handling
%include "exception.i"

// Director exception handling
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Teuchos::RCPs typemaps
%teuchos_rcp(NOX::Epetra::Interface::Required)
%teuchos_rcp(NOX::Epetra::Interface::Jacobian)
%teuchos_rcp(NOX::Epetra::Interface::Preconditioner)

///////////////////////
// NOX_Utils support //
///////////////////////
%import "NOX_Utils.i"

///////////////////////////////////////////
// NOX_Epetra_Interface_Required support //
///////////////////////////////////////////
%feature("autodoc",
"computeF(self, Epetra.Vector x, Epetra.Vector F, FillType flag) -> bool

  Virtual method in C++ that is intended to be overridden by user.
  This method defines the nonlinear function to be solved.  Arguments
  x and F will be provided as numpy-hybrid Epetra.Vector objects.
  Return True if the computation is successful.

  It is strongly recommended that computeF() not raise any exceptions,
  accidental or otherwise.  This can be prevented by wrapping your
  algorithm in a try block:

    try:
      # Your code here...
    except Exception, e:
      print 'Python exception raised in computeF():'
      print e
      return False

  By returning False, you tell NOX that computeF() was unsuccessful.
")
NOX::Epetra::Interface::Required::computeF;
%feature("director") NOX::Epetra::Interface::Required;
%include "NOX_Epetra_Interface_Required.H"

///////////////////////////////////////////
// NOX_Epetra_Interface_Jacobian support //
///////////////////////////////////////////
%feature("director") NOX::Epetra::Interface::Jacobian;
%include "NOX_Epetra_Interface_Jacobian.H"

/////////////////////////////////////////////////
// NOX_Epetra_Interface_Preconditioner support //
/////////////////////////////////////////////////
%feature("director") NOX::Epetra::Interface::Preconditioner;
%include "NOX_Epetra_Interface_Preconditioner.H"
