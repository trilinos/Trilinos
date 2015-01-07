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

////////////////////////////////////////////////////////
// Typemaps for NOX::Petsc::Interface virtual methods //
////////////////////////////////////////////////////////

%typemap(directorin) Vec &
{
  $input = PyPetscVec_New($1);
}

%typemap(directorin) const Vec &
{
  $input = PyPetscVec_New($1);
}

%typemap(directorin) Mat &
{
  $input = PyPetscMat_New($1);
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
%teuchos_rcp(NOX::Petsc::Interface)
%feature("director") NOX::Petsc::Interface;
%include "NOX_Petsc_Interface.H"
