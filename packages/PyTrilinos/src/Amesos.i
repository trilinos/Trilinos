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
"The Amesos module allows access to The Trilinos package Amesos.  Note
that the 'Amesos_' prefix has been stripped from all Amesos objects,
but that if imported with 'from PyTrilinos import Amesos', these
objects exist in the 'Amesos' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.

The most important classes of the Amesos module are:
- The factory class, Amesos.Factory
- The LAPACK interface, Amesos.Lapack
- The KLU interace, Amesos.Klu
- The UMFPACK interace, Amesos.Umfpack
- The SCALAPACK interace, Amesos.Scalapack
- The SuperLU interace, Amesos.Superlu
- The SuperLU_DIST interace, Amesos.Superludist
- The DSCPACK interace, Amesos.Dscpack
- The MUMPS interace, Amesos.Mumps

Each specific interface may require Amesos to be configured with the
appropriate --enable-<interface>. (Note that LAPACK and KLU are enabled
by default.)

For examples of usage, please consult the python/example subdirectory.
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
#include "Amesos_ConfigDefs.h"

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_DataAccess.h"
#include "Epetra_FEVbrMatrix.h"

// Teuchos includes
#include "Teuchos_RefCountPtr.hpp"

// Teuchos Python utility code
#include "Teuchos_PythonParameter.h"

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
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

%}

// External Trilinos packages
%import "Teuchos.i"
%import "Epetra.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Rename directives for Amesos
%rename(BaseSolver ) Amesos_BaseSolver;
%rename(Factory    ) Amesos;
%rename(Klu        ) Amesos_Klu;
%rename(Lapack     ) Amesos_Lapack;
%rename(Umfpack    ) Amesos_Umfpack;
%rename(Scalapack  ) Amesos_Scalapack;
%rename(Taucs      ) Amesos_Taucs;
%rename(Pardiso    ) Amesos_Pardiso;
%rename(Superlu    ) Amesos_Superlu;
%rename(Superludist) Amesos_Superludist;
%rename(Mumps      ) Amesos_Mumps;
%rename(Dscpack    ) Amesos_Dscpack;

// SWIG library includes
%include "stl.i"

/////////////////////////
// Amesos core support //
/////////////////////////
%include "Amesos_config.h"
%include "Amesos_ConfigDefs.h"
%include "Amesos.h"
%include "Amesos_BaseSolver.h"
%extend Amesos_BaseSolver 
{
  string __str__() {
    stringstream os;
    os << "*** Amesos_BaseSolver ***";
    return os.str();
  }
  void __del__() {
    delete self;
  }
}

///////////////////////////
// Amesos LAPACK support //
///////////////////////////
#ifdef HAVE_AMESOS_LAPACK
%include "Amesos_Lapack.h"
#endif

////////////////////////
// Amesos KLU support //
////////////////////////
#ifdef HAVE_AMESOS_KLU
%include "Amesos_Klu.h"
#endif

////////////////////////////
// Amesos UMFPACK support //
////////////////////////////
#ifdef HAVE_AMESOS_UMFPACK
%include "Amesos_Umfpack.h"
#endif

//////////////////////////////
// Amesos ScaLAPACK support //
//////////////////////////////
#ifdef HAVE_AMESOS_SCALAPACK
%include "Amesos_Scalapack.h"
#endif

//////////////////////////
// Amesos Taucs support //
//////////////////////////
#ifdef HAVE_AMESOS_TAUCS
%include "Amesos_Taucs.h"
#endif

////////////////////////////
// Amesos Pardiso support //
////////////////////////////
#ifdef HAVE_AMESOS_PARDISO
%include "Amesos_Pardiso.h"
#endif

////////////////////////////
// Amesos SuperLU support //
////////////////////////////
#ifdef HAVE_AMESOS_SUPERLU
%include "Amesos_Superlu.h"
#endif

////////////////////////////////
// Amesos SuperLUDist support //
////////////////////////////////
#ifdef HAVE_AMESOS_SUPERLUDIST
%include "Amesos_Superludist.h"
#endif

//////////////////////////
// Amesos MUMPS support //
//////////////////////////
#ifdef HAVE_AMESOS_MUMPS
%include "Amesos_Mumps.h"
#endif
