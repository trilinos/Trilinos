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

%define %amesos_docstring
"
PyTrilinos.Amesos is the python interface to the Trilinos direct
linear solver package Amesos:

    http://trilinos.sandia.gov/packages/amesos

The purpose of Amesos is to provide a common interface to a variety of
third-party direct solvers, made compatible with PyTrilinos.Epetra.
Note that the C++ version of Amesos uses the prefix 'Amesos_', which
has been stripped from the python implementation.

The most important classes of the Amesos module are:

    * Factory      - Factory class
    * Lapack       - LAPACK interface
    * Klu          - KLU interface
    * Umfpack      - UMFPACK interface
    * Scalapack    - SCALAPACK interface
    * Superlu      - SuperLU interface
    * Superludist  - SuperLU_DIST interface
    * Dscpack      - DSCPACK interface
    * Mumps        - MUMPS interface

Use dir(Amesos) to see what specific interfaces have been enabled on
your platform.  For examples of usage, please consult the examples
subdirectory of the PyTrilinos package, scripts exAmesos_Simple.py and
exAmesos_Factory.py.
"
%enddef

%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %amesos_docstring) Amesos

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
#include "Amesos_ConfigDefs.h"

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_Object.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_DataAccess.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_SerialDistributor.h"
#endif

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_RefCountPtrDecl.hpp"
#include "PyTrilinos_Teuchos_Util.h"
#endif

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#ifdef HAVE_EPETRA
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#endif

%}

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include Amesos documentation
%include "Amesos_dox.i"

// SWIG library includes
%include "stl.i"

// External Trilinos packages
#ifdef HAVE_TEUCHOS
%import "Teuchos.i"
#endif
#ifdef HAVE_EPETRA
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
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

//////////////////////////
// Amesos configuration //
//////////////////////////
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
%include "Amesos_config.h"
%include "Amesos_ConfigDefs.h"

//////////////////////////////
// Amesos (Factory) support //
//////////////////////////////
%rename(Factory) Amesos;
%newobject Amesos::Create;
%include "Amesos.h"

///////////////////////////////
// Amesos BaseSolver support //
///////////////////////////////
%rename(BaseSolver) Amesos_BaseSolver;
%include "Amesos_BaseSolver.h"

///////////////////////////
// Amesos LAPACK support //
///////////////////////////
#ifdef HAVE_AMESOS_LAPACK
%rename(Lapack) Amesos_Lapack;
%include "Amesos_Lapack.h"
#endif

////////////////////////
// Amesos KLU support //
////////////////////////
#ifdef HAVE_AMESOS_KLU
%rename(Klu) Amesos_Klu;
%include "Amesos_Klu.h"
#endif

////////////////////////////
// Amesos UMFPACK support //
////////////////////////////
#ifdef HAVE_AMESOS_UMFPACK
%rename(Umfpack) Amesos_Umfpack;
%include "Amesos_Umfpack.h"
#endif

//////////////////////////////
// Amesos ScaLAPACK support //
//////////////////////////////
#ifdef HAVE_AMESOS_SCALAPACK
%rename(Scalapack) Amesos_Scalapack;
%include "Amesos_Scalapack.h"
#endif

//////////////////////////
// Amesos Taucs support //
//////////////////////////
#ifdef HAVE_AMESOS_TAUCS
%rename(Taucs) Amesos_Taucs;
%include "Amesos_Taucs.h"
#endif

////////////////////////////
// Amesos Pardiso support //
////////////////////////////
#ifdef HAVE_AMESOS_PARDISO
%rename(Pardiso) Amesos_Pardiso;
%include "Amesos_Pardiso.h"
#endif

////////////////////////////
// Amesos SuperLU support //
////////////////////////////
#ifdef HAVE_AMESOS_SUPERLU
%rename(Superlu) Amesos_Superlu;
%include "Amesos_Superlu.h"
#endif

////////////////////////////////
// Amesos SuperLUDist support //
////////////////////////////////
#ifdef HAVE_AMESOS_SUPERLUDIST
%rename(Superludist) Amesos_Superludist;
%include "Amesos_Superludist.h"
#endif

//////////////////////////
// Amesos MUMPS support //
//////////////////////////
#ifdef HAVE_AMESOS_MUMPS
%rename(Mumps) Amesos_Mumps;
%include "Amesos_Mumps.h"
#endif

////////////////////////////
// Amesos DSCPACK support //
////////////////////////////
#ifdef HAVE_AMESOS_DSCPACK
%rename(Dscpack) Amesos_Dscpack;
%include "Amesos_Dscpack.h"
#endif

// Turn off the exception handling
%exception;
