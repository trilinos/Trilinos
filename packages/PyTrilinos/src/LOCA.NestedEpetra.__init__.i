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

%define %loca_epetra_docstring
"
PyTrilinos.LOCA.Epetra is the python interface to namespace Epetra of
the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Epetra is to provide ***.  The python version of
LOCA.Epetra supports the following classes:

    * Group  -

Any other notes about the package as a whole. . . .
"
%enddef

%module(package      = "PyTrilinos.LOCA.NestedEpetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_epetra_docstring) __init__

%{
// System includes
#include <vector>

// PyTrilinos includes
#include "PyTrilinos_Teuchos_Util.h"
#include "PyTrilinos_Epetra_Util.h"

// Local Epetra includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"

// Epetra includes
#include "Epetra_DLLExportMacro.h"
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

// NOX includes
#include "NOX.H"
#include "NOX_Epetra_Group.H"

// LOCA includes
#include "LOCA.H"
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"
#include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_Constraint.H"
#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
#include "LOCA_Epetra.H"
#include "LOCA_Epetra_Group.H"

// Namespace flattening
using Teuchos::RCP;
using Teuchos::rcp;
%}

%ignore *::operator=;

// SWIG library includes
%include "stl.i"

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

%include "Epetra_DLLExportMacro.h"

// Trilinos interface support
%import "Teuchos.i"
%include "Epetra_Base.i"   // For PyExc_EpetraError
%import "Epetra.i"
%import "NOX.Abstract.i"
%import "NOX.NestedEpetra.__init__.i"
%import "NOX.NestedEpetra.Interface.i"
%import "LOCA.__init__.i"
%import "LOCA.MultiContinuation.i"
%import "LOCA.Homotopy.i"
%import "LOCA.Hopf.MinimallyAugmented.i"
%import "LOCA.Pitchfork.MinimallyAugmented.i"
%import "LOCA.TurningPoint.MinimallyAugmented.i"
%import "LOCA.Abstract.i"
%import "LOCA.NestedEpetra.Interface.i"

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

/////////////////////////
// LOCA Epetra support //
/////////////////////////

#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
%include "LOCA_Epetra.H"

//////////////////////////////
// LOCA.Epetra.Group support //
//////////////////////////////

// temporarily ignore conflict-causing constructor.  TODO: fix this issue
%ignore LOCA::Epetra::Group::Group(Teuchos::RCP< LOCA::GlobalData > const &,Teuchos::ParameterList &,Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree > const &,NOX::Epetra::Vector &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,Teuchos::RCP< NOX::Epetra::LinearSystem > const &,LOCA::ParameterVector const &);

%teuchos_rcp(LOCA::Epetra::Group)
%include "LOCA_Epetra_Group.H"
