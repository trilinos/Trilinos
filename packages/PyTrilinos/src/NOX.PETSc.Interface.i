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

%define %nox_petsc_interface_docstring
"
PyTrilinos.NOX.PETSc.Interface is the python interface to the
PETSc::Interface namespace of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.PETSc.Interface is to provide base classes the
user should derive from in order to define the nonlinear function to
be solved, and if needed, its Jacobian and the desired preconditioner.

NOX.Epetra.Interface provides the following user-level class:

    * Interface       - Required class for computing the nonlinear function
"
%enddef

%module(package      = "PyTrilinos.NOX.PETSc",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_petsc_interface_docstring) Interface

%{
// PyTrilinos includes
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.hpp"

// NumPy includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
// #include "PyTrilinos_Teuchos_Util.hpp"

// PETSc include
#include "petsc.h"

// NOX::Petsc::Interface includes
#include "NOX_Petsc_Interface.H"
%}

// Include the PETSc4Py SWIG interface file
%include "petsc4py/petsc4py.i"

// Include NOX documentation
%include "NOX_dox.i"

// STL support
%include "stl.i"

// Trilinos module imports
%import "Teuchos.i"

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

/////////////////////////////////
// NOX_Petsc_Interface support //
/////////////////////////////////
// %feature("autodoc",
// "computeF(self, Epetra.Vector x, Epetra.Vector F, FillType flag) -> bool

//   Virtual method in C++ that is intended to be overridden by user.
//   This method defines the nonlinear function to be solved.  Arguments
//   x and F will be provided as numpy-hybrid Epetra.Vector objects.
//   Return True if the computation is successful.

//   It is strongly recommended that computeF() not raise any exceptions,
//   accidental or otherwise.  This can be prevented by wrapping your
//   algorithm in a try block:

//     try:
//       # Your code here...
//     except Exception, e:
//       print 'Python exception raised in computeF():'
//       print e
//       return False

//   By returning False, you tell NOX that computeF() was unsuccessful.
// ")
// NOX::Epetra::Interface::Required::computeF;
%teuchos_rcp(NOX::Petsc::Interface::Interface)
%feature("director") NOX::Petsc::Interface::Interface;
%include "NOX_Petsc_Interface.H"
