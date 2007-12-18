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

// The swig interface file Teucos.i defines typemaps that greatly
// expand the usefulness of the RCP class.  The Epetra.i swig
// interface file defines new classes that inherit from Epetra
// array-like classes (such as Epetra_MultiVector or Epetra_Vector)
// that also behave like NumPy arrays, making them much more usable in
// python.

// The purpose of this swig interface file is provide for the cases
// that involve an RCP<...> wrapped around an Epetra_Vector or other
// array-like class.  In this case, what the python user will expect
// is an enhanced Epetra array-like object. To use this interface
// file, use the directive
//
// %include "Teuchos_Epetra.i"
//
// It is not necessary (nor harmful) to also use '%import "Teuchos.i"'
// and/or '%import "Epetra.i"'.

// Forward declarations
%{
// PyObject * PyExc_EpetraError;
// std::string Epetra_Object___str__(Epetra_Object*);
// void        Epetra_Object_Print(  Epetra_Object*,PyObject*pf=NULL);

// Namespace flattening
using Teuchos::RCP;
%}

// Imports
%import "Teuchos.i"
%import "Epetra.i"

// Define the macro that defines the typemap
%define %teuchos_rcp_epetra_array_typemaps(ClassName)
%typemap(out) RCP<Epetra_##ClassName>
{
  if (Teuchos::is_null($1)) $result = Py_BuildValue("");
  else
  {
    Epetra_NumPy##ClassName * npa = new Epetra_NumPy##ClassName(*$1);
    $result = SWIG_NewPointerObj(npa, $descriptor(Epetra_NumPy##ClassName*), 1);
  }
}
%typemap(varout) RCP<Epetra_##ClassName>
{
  if (Teuchos::is_null($1)) $result = Py_BuildValue("");
  else
  {
    Epetra_NumPy##ClassName * npa = new Epetra_NumPy##ClassName(*$1);
    $result = SWIG_NewPointerObj(npa, $descriptor(Epetra_NumPy##ClassName*), 1);
  }
}

// %ignore RCP< Epetra_##ClassName >::access_node() const;
// %template (RCP_Epetra_##ClassName)
//   RCP< Epetra_##ClassName >;
%enddef

// Implement the macro for concrete classes
%teuchos_rcp_epetra_array_typemaps(IntVector           )
%teuchos_rcp_epetra_array_typemaps(MultiVector         )
%teuchos_rcp_epetra_array_typemaps(Vector              )
%teuchos_rcp_epetra_array_typemaps(FEVector            )
%teuchos_rcp_epetra_array_typemaps(IntSerialDenseMatrix)
%teuchos_rcp_epetra_array_typemaps(IntSerialDenseVector)
%teuchos_rcp_epetra_array_typemaps(SerialDenseMatrix   )
%teuchos_rcp_epetra_array_typemaps(SerialDenseVector   )
