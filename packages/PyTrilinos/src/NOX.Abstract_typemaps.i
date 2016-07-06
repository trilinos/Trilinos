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

%{
// Teuchos includes
#include "Teuchos_RCP.hpp"

// PyTrilinos includes
#include "PyTrilinos_config.h"

// Epetra includes
#ifdef HAVE_NOX_EPETRA
#include "Epetra_Vector.h"
#endif

// NOX includes
#include "NOX_Abstract_Vector.H"
#ifdef HAVE_NOX_PETSC
#include "NOX_Petsc_Vector.H"
#endif
%}

#ifdef HAVE_NOX_PETSC
%include "petsc4py/petsc4py.i"
#endif

/////////////////////////////////////////////////////////////////////
// *** Utility functions for downcasting NOX::Abstract classes *** //
/////////////////////////////////////////////////////////////////////

%fragment("NOX_Abstract_downcast", "header")
{
  /////////////////////////////////////////////////////////////////////
  // A utility function that attempts to downcast a
  // NOX::Abstract::Vector to a viable subclass and return as a wrapped
  // Python object.  Viable subclasses currently include
  // NOX::Epetra::Vector and NOX::Petsc::Vector.
  /////////////////////////////////////////////////////////////////////

  PyObject *
    downcast_NOX_Abstract_Vector(NOX::Abstract::Vector & nav,
                                 bool owner)
  {
    static swig_type_info * swig_NAV_ptr =
      SWIG_TypeQuery("NOX::Abstract::Vector *");

#ifdef HAVE_NOX_EPETRA
    // Try to downcast to a NOX::Epetra::Vector
    static swig_type_info * swig_EV_ptr =
      SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
    NOX::Epetra::Vector * nevResult =
      dynamic_cast< NOX::Epetra::Vector * >(&nav);
    if (nevResult != NULL)
    {
      Teuchos::RCP< Epetra_Vector > *smartresult = new
        Teuchos::RCP< Epetra_Vector >(new
          Epetra_Vector(View,
                        nevResult->getEpetraVector(),
                        0),
                                      owner);
      return SWIG_NewPointerObj(SWIG_as_voidptr(smartresult),
                                swig_EV_ptr,
                                SWIG_POINTER_OWN);
    }
#endif

#ifdef HAVE_NOX_PETSC
    // Try to downcast to a NOX::Petsc::Vector
    NOX::Petsc::Vector * npvResult = dynamic_cast< NOX::Petsc::Vector * >(&nav);
    if (npvResult != NULL)
    {
      return PyPetscVec_New(npvResult->getPetscVector());
    }
#endif

    // No downcasts worked, so return as a Python wrapped
    // NOX::Abstract::Vector
    return SWIG_NewPointerObj(SWIG_as_voidptr(&nav),
                              swig_NAV_ptr,
                              SWIG_POINTER_OWN);
  }

  /////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////
  // A utility function that attempts to downcast a const
  // NOX::Abstract::Vector to a viable subclass and return as a wrapped
  // Python object.  Viable subclasses currently include const
  // NOX::Epetra::Vector and const NOX::Petsc::Vector.
  /////////////////////////////////////////////////////////////////////

  PyObject *
    downcast_NOX_Abstract_Vector(const NOX::Abstract::Vector & nav,
                                 bool owner)
  {
    static swig_type_info * swig_NAV_ptr =
      SWIG_TypeQuery("const NOX::Abstract::Vector *");

#ifdef HAVE_NOX_EPETRA
    // Try to downcast to a NOX::Epetra::Vector
    static swig_type_info * swig_EV_ptr =
      SWIG_TypeQuery("Teuchos::RCP< const Epetra_Vector > *");
    const NOX::Epetra::Vector * nevResult =
      dynamic_cast< const NOX::Epetra::Vector * >(&nav);
    if (nevResult != NULL)
    {
      Teuchos::RCP< const Epetra_Vector > * smartresult = new
        Teuchos::RCP< const Epetra_Vector >(new
          const Epetra_Vector(View,
                              nevResult->getEpetraVector(),
                              0),
                                            owner);
      return SWIG_NewPointerObj(SWIG_as_voidptr(smartresult),
                                swig_EV_ptr,
                                SWIG_POINTER_OWN);
    }
#endif

#ifdef HAVE_NOX_PETSC
    // Try to downcast to a const NOX::Petsc::Vector
    const NOX::Petsc::Vector * npvResult =
      dynamic_cast< const NOX::Petsc::Vector * >(&nav);
    if (npvResult != NULL)
    {
      return PyPetscVec_New(npvResult->getPetscVector());
    }
#endif

    // No downcasts worked, so return as a Python wrapped
    // NOX::Abstract::Vector
    return SWIG_NewPointerObj(SWIG_as_voidptr(&nav),
                              swig_NAV_ptr,
                              SWIG_POINTER_OWN);
  }
}

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// *** Typemaps ***
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////

// Take an output NOX::Abstract::Vector and downcast to either
// Epetra.Vector or PETSc.Vec, if possible, and convert to Python
%typemap(out, fragment="NOX_Abstract_downcast")
  NOX::Abstract::Vector &
{
  $result = downcast_NOX_Abstract_Vector(*$1, bool($owner));
}

/////////////////////////////////////////////////////////////////////

#ifdef HAVE_NOX_EPETRA

/////////////////////////////////////////////////////////////////////

// Make Epetra_Vector and NOX::Epetra::Vector input arguments
// interchangeable
%typemap(in) NOX::Epetra::Vector &
(void* argp=0,
 int res=0,
 Teuchos::RCP< Epetra_Vector > tempshared,
 bool cleanup=false)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res))
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input,
				&argp,
				$descriptor(Teuchos::RCP< Epetra_Vector > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp)
    {
      %argument_nullref("$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< Epetra_Vector > *);
      delete %reinterpret_cast(argp, Teuchos::RCP< Epetra_Vector > *);
      $1 = new NOX::Epetra::Vector(Teuchos::rcp_dynamic_cast< Epetra_Vector >(tempshared),
				   NOX::Epetra::Vector::CreateView);
      cleanup = true;
    }
    else
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< Epetra_Vector > *);
      $1 = new NOX::Epetra::Vector(Teuchos::rcp_dynamic_cast< Epetra_Vector >(tempshared),
				   NOX::Epetra::Vector::CreateView);
      cleanup = true;
    }
  }
  else
  {
    $1 = %reinterpret_cast(argp, NOX::Epetra::Vector*);
  }
}

/////////////////////////////////////////////////////////////////////

%typecheck(1190) NOX::Epetra::Vector &
{
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, 0, $descriptor, 0)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtrAndOwn($input, 0,
			       $descriptor(Teuchos::RCP< Epetra_Vector > *),
			       %convertptr_flags, 0)) ? 1 : 0;
}

/////////////////////////////////////////////////////////////////////

%typecheck(1190) const NOX::Epetra::Vector &
{
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, 0, $descriptor, 0)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtrAndOwn($input, 0,
			       $descriptor(Teuchos::RCP< Epetra_Vector > *),
			       %convertptr_flags, 0)) ? 1 : 0;
}

/////////////////////////////////////////////////////////////////////

%typemap(freearg) NOX::Epetra::Vector &
{
  if (cleanup$argnum) delete $1;
}

/////////////////////////////////////////////////////////////////////

// Convert NOX::Epetra::LinearSystem objects to
// NOX::Epetra::LinearSystemAztecOO
%typemap(out) Teuchos::RCP< NOX::Epetra::LinearSystem >
(NOX::Epetra::LinearSystem*        nelsPtr     = NULL,
 NOX::Epetra::LinearSystemAztecOO* nelsaResult = NULL)
{
  nelsPtr = $1.get();
  nelsaResult = dynamic_cast< NOX::Epetra::LinearSystemAztecOO*>(nelsPtr);
  if (nelsaResult == NULL)
  {
    //If we cannot downcast then return the NOX::Epetra::LinearSystem
    %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystem > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *smartresult =
      new Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO >(nelsaResult);
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *),
				   SWIG_POINTER_OWN));
  }
}

/////////////////////////////////////////////////////////////////////

%typemap(out) Teuchos::RCP< const NOX::Epetra::LinearSystem >
(const NOX::Epetra::LinearSystem*        nelsPtr     = NULL,
 const NOX::Epetra::LinearSystemAztecOO* nelsaResult = NULL)
{
  nelsPtr = $1.get();
  nelsaResult = dynamic_cast< const NOX::Epetra::LinearSystemAztecOO*>(nelsPtr);
  if (nelsaResult == NULL)
  {
    //If we cannot downcast then return the NOX::Epetra::LinearSystem
    %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystem > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    Teuchos::RCP< const NOX::Epetra::LinearSystemAztecOO > *smartresult =
      new Teuchos::RCP< const NOX::Epetra::LinearSystemAztecOO >(nelsaResult);
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Epetra::LinearSystemAztecOO > *),
				   SWIG_POINTER_OWN));
  }
}

#endif
