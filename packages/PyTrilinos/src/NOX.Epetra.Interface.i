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
#include "numpy_include.hpp"

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

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
#include "PyTrilinos_Epetra_Util.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"

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
  catch(PyTrilinos::PythonException & e)
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
