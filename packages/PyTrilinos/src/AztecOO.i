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

%define %aztecoo_docstring
"
PyTrilinos.AztecOO is the python interface to the Trilinos iterative
linear solver package AztecOO:

    http://trilinos.sandia.gov/packages/aztecoo

AztecOO is the object-oriented interface to Aztec, Sandia's venerable
Krylov-space linear system solver package.  Note that the C++ version
of AztecOO uses the prefix 'AztecOO_' which has been stripped from the
python version.  AztecOO requires the Epetra module. The IFPACK and ML
modules can extend the preconditioning capabilities of AztecOO.

AztecOO has a single class:

    * AztecOO  - Object-oriented interface to Aztec package

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exAztecOO.py
    * exAztecOO_Operator.py
    * exAztecOO_RowMatrix.py
    * exAztecOO_BasicRowMatrix.py
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %aztecoo_docstring) AztecOO

#if SWIG_VERSION >= 0x030000
%feature("flatnested");
#else
// Handle the AztecOO::MatrixData and AztecOO::OperatorData nested
// structs by defining them exclusively for SWIG as though they were
// not nested
struct MatrixData
{
  Epetra_RowMatrix * A;
  Epetra_Vector * X;
  Epetra_Vector * Y;
  Epetra_Vector * SourceVec;
  Epetra_Vector * TargetVec;
  MatrixData(Epetra_RowMatrix * inA = 0, Epetra_Vector * inX = 0,
	     Epetra_Vector * inY = 0, Epetra_Vector * inSourceVec = 0,
	     Epetra_Vector * inTargetVec = 0);
  ~MatrixData();
};
%nestedworkaround AztecOO::MatrixData;

struct OperatorData
{
  Epetra_Operator * A;
  Epetra_Vector * X;
  Epetra_Vector * Y;
  OperatorData(Epetra_Operator * inA = 0, Epetra_Vector * inX = 0,
	       Epetra_Vector * inY = 0);
  ~OperatorData();
};
%nestedworkaround AztecOO::OperatorData;
#endif

%{
// System include files
#include <iostream>
#include <sstream>
#include <vector>

// Configuration include files
#include "PyTrilinos_config.h"
#include "AztecOO_ConfigDefs.h"

// Optional Teuchos support
#ifdef HAVE_AZTECOO_TEUCHOS
#include "PyTrilinos_Teuchos_Headers.hpp"
#endif

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"

// NumPy include
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#endif

// AztecOO include files
#include "PyTrilinos_AztecOO_Headers.hpp"
%}

// Auto-documentation feature
%feature("autodoc", "1");

%include "PyTrilinos_config.h"

// AztecOO enumerated types support
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
%include "AztecOO_config.h"
%include "az_aztec_defs.h"

// Include AztecOO documentation
%include "AztecOO_dox.i"

// Include the NumPy typemaps
%include "numpy.i"

// Include the standard exception handlers
%include "exception.i"

// External Trilinos interface imports
#ifdef HAVE_PYTRILINOS_EPETRA
%import "Epetra.i"
#endif
#ifdef HAVE_AZTECOO_TEUCHOS
%import "Teuchos.i"
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
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(int errCode)
  {
    PyErr_Format(PyExc_RuntimeError, "Error code = %d\nSee stderr for details",
		 errCode);
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// Macro for methods that return C arrays
%define %aztecoo_return_array(className,methodName,type,typeName,length)
%extend className
{
  PyObject * methodName() const
  {
    npy_intp dims[ ] = { (npy_intp) length };
    return PyArray_SimpleNewFromData(1, dims, typeName, (void*)self->methodName());
  }
}
%ignore className::methodName() const;
%enddef

/////////////////////////////
// AztecOO Version support //
/////////////////////////////
%include "AztecOO_Version.h"
%pythoncode
%{
__version__ = AztecOO_Version().split()[3]
%}

/////////////////////
// AztecOO support //
/////////////////////
%ignore AztecOO::GetAllAztecStatus(double*);
%aztecoo_return_array(AztecOO, GetAllAztecOptions, int,    NPY_INT,    AZ_OPTIONS_SIZE)
%aztecoo_return_array(AztecOO, GetAllAztecParams,  double, NPY_DOUBLE, AZ_PARAMS_SIZE )
%aztecoo_return_array(AztecOO, GetAztecStatus,     double, NPY_DOUBLE, AZ_STATUS_SIZE )
%include "AztecOO.h"
// We've fooled SWIG into thinking that MatrixData and OperatorData
// are global structs, so now we need to trick the C++ compiler into
// understanding these apparent global types.
%{
typedef AztecOO::MatrixData   MatrixData;
typedef AztecOO::OperatorData OperatorData;
%}

// Turn off the exception handling
%exception;
