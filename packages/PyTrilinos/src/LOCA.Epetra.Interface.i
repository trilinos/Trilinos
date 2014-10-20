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

%define %loca_epetra_interface_docstring
"
PyTrilinos.LOCA.Epetra.Interface is the python interface to namespace
Epetra::Interface of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra.Interface is to provide a concrete Epetra
implementation of LOCA interfaces.  The python version of
LOCA.Epetra.Interface supports the following classes:

    * Required                 - Provides a set of interfaces for users to
                                 provide information about the nonlinear
                                 problem to LOCA
    * MassMatrix               - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for the coefficients of
                                 time dependent terms
    * TimeDependent            - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for computing the shifted
                                 matrix
    * TimeDependentMatrixFree  - Used by LOCA.Epetra.Group to provide a link
                                 to the external code for applying the shifted
                                 matrix in a matrix-free setting
"
%enddef

%module(package      = "PyTrilinos.LOCA.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_interface_docstring) Interface

%{
// NumPy includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// PyTrilinos includes
#include "PyTrilinos_Teuchos_Util.hpp"
#include "PyTrilinos_Epetra_Util.hpp"

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// Local Epetra includes
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"

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

// NOX include
#include "NOX_Epetra_Interface_Required.H"

// LOCA::Epetra::Interface includes
#include "LOCA_Epetra_Interface_Required.H"
#include "LOCA_Epetra_Interface_MassMatrix.H"
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"
%}

// General ignore directives
%ignore *::operator=;

// STL support
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

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
  catch(PyTrilinos::PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  catch(int errCode)
  {
    PyErr_Format(PyExc_EpetraError, "Error code = %d\nSee stderr for details", errCode);
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
%teuchos_rcp(LOCA::Epetra::Interface::Required)
%teuchos_rcp(LOCA::Epetra::Interface::MassMatrix)
%teuchos_rcp(LOCA::Epetra::Interface::TimeDependent)
%teuchos_rcp(LOCA::Epetra::Interface::TimeDependentMatrixFree)

///////////////////////
// NOX_Utils support //
///////////////////////
%import "NOX_Utils.i"


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
LOCA::Epetra::Interface::Required::computeF;

%import "NOX.Epetra.Interface.i"

%feature("director") LOCA::Epetra::Interface::Required;
%include "LOCA_Epetra_Interface_Required.H"

%feature("director") LOCA::Epetra::Interface::MassMatrix;
%include "LOCA_Epetra_Interface_MassMatrix.H"

%feature("director") LOCA::Epetra::Interface::TimeDependent;
%include "LOCA_Epetra_Interface_TimeDependent.H"

%feature("director") LOCA::Epetra::Interface::TimeDependentMatrixFree;
// The following #define is to change the name of LOCA method
// arguments that conflict with a SWIG director method argument
#define result loca_result
%include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"
#undef result
