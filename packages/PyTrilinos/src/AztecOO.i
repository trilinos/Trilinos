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

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_SYS_TIME_H
#undef HAVE_SYS_TIME_H
#endif
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "AztecOO_ConfigDefs.h"

// Optional Teuchos support
#ifdef HAVE_AZTECOO_TEUCHOS
#include "Teuchos_PythonParameter.h"
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
#include "Epetra_FEVector.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"
#include "Epetra_SerialDistributor.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

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

// AztecOO includes
#include "AztecOO.h"
#include "AztecOO_Version.h"

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
#ifdef HAVE_EPETRA
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
