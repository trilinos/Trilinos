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
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"
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
#include "numpy_include.hpp"
#ifdef HAVE_EPETRA
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyMultiVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
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
