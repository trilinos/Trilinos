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


#ifdef HAVE_NOX_EPETRA
// Downcast NOX::Abstract::Vector return arguments to Epetra.Vectors,
// if possible
%typemap(out) NOX::Abstract::Vector &
(NOX::Epetra::Vector* nevResult = NULL)
{
  nevResult = dynamic_cast<NOX::Epetra::Vector*>($1);
  if (nevResult == NULL)
  {
    // If we cannot downcast, then return the NOX::Abstract::Vector
    $result = SWIG_NewPointerObj((void*)&$1, $descriptor, SWIG_POINTER_OWN);
  }
  else
  {
    Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *smartresult = new
      Teuchos::RCP< PyTrilinos::Epetra_NumPyVector >(new PyTrilinos::Epetra_NumPyVector(View,
								nevResult->getEpetraVector(),
								0), bool($owner));
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *),
				   SWIG_POINTER_OWN));
  }
}

%typemap(out) Teuchos::RCP< const NOX::Abstract::Vector >
(Teuchos::RCP< const NOX::Epetra::Vector > nevResult)
{
  nevResult = Teuchos::rcp_dynamic_cast< const NOX::Epetra::Vector >(*(&$1));
  if (nevResult.is_null())
  {
    // If we cannot downcast, then return the NOX::Abstract::Vector
    Teuchos::RCP< const NOX::Abstract::Vector > *smartresult = $1.is_null() ? 0 :
      new Teuchos::RCP< const NOX::Abstract::Vector >($1);
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< NOX::Abstract::Vector > *),
				   SWIG_POINTER_OWN));
  }
  else
  {
    Teuchos::RCP< const PyTrilinos::Epetra_NumPyVector > *smartresult = new
      Teuchos::RCP< const PyTrilinos::Epetra_NumPyVector >(new PyTrilinos::Epetra_NumPyVector(View,
								      (*nevResult).getEpetraVector(),
								      0), bool($owner));
    %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
				   $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *),
				   SWIG_POINTER_OWN));
  }
}


// Make Epetra_Vector and NOX::Epetra::Vector input arguments
// interchangeable
%typemap(in) NOX::Epetra::Vector &
(void* argp=0, int res=0, Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > tempshared,
 bool cleanup=false)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res))
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input,
				&argp,
				$descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *),
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
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *);
      delete %reinterpret_cast(argp, Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *);
      $1 = new NOX::Epetra::Vector(Teuchos::rcp_dynamic_cast< Epetra_Vector >(tempshared),
				   NOX::Epetra::Vector::CreateView);
      cleanup = true;
    }
    else
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *);
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
%typecheck(1190) NOX::Epetra::Vector &
{
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, 0, $descriptor, 0)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtrAndOwn($input, 0,
			       $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPyVector > *),
			       %convertptr_flags, 0)) ? 1 : 0;
}
%typemap(freearg) NOX::Epetra::Vector &
{
  if (cleanup$argnum) delete $1;
}

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
