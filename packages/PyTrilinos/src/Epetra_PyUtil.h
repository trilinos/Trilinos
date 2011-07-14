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

#ifndef EPETRA_PYUTIL_H
#define EPETRA_PYUTIL_H

// Include the Python prototypes
#include "Python.h"

// Python developers improved the const-correctness of char* variables
// in the C-API with the advent of version 2.5.
#if PY_VERSION_HEX >= 0x02050000
#define CONST const
#else
#define CONST
#endif

// PyTrilinos include
#include "PyTrilinos_config.h"

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"
#endif

// Epetra includes
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"


////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////

// Given a pointer to an Epetra_MultiVector, convert to a python
// object and return the pointer.  Attempt to downcast to an
// Epetra_NumPyMultiVector.
#ifdef HAVE_TEUCHOS
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< Epetra_MultiVector > *emv);
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< const Epetra_MultiVector > *emv);
#else
PyObject *
convertEpetraMultiVectorToPython(const Epetra_MultiVector * emv);
#endif

// Given a pointer to an Epetra_Vector, convert to a python object and
// return the pointer.  Attempt to downcast to an Epetra_NumPyVector.
#ifdef HAVE_TEUCHOS
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< Epetra_Vector > *ev);
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< const Epetra_Vector > *ev);
#else
PyObject *
convertEpetraVectorToPython(const Epetra_Vector * ev);
#endif

// Given a pointer to an Epetra_Operator, convert to a python
// object and return the pointer.  Attempt to downcast to any one of
// the many Epetra classes that derive from Epetra_Operator.
#ifdef HAVE_TEUCHOS
PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< Epetra_Operator > *eo);
PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< const Epetra_Operator > *eo);
#else
PyObject *
convertEpetraOperatorToPython(const Epetra_Operator * eo, int cnvt_flags=0);
#endif

////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

// Given a const Epetra_BlockMap, return a reference counted pointer
// to a const Epetra_Map.  If the downcast cannot be performed, throw a
// PythonException.
Teuchos::RCP<const Epetra_Map>
getEpetraMapPtrFromEpetraBlockMap(const Epetra_BlockMap & ebm);

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_Vector value of the attribute.  If
// the attribute does not exist or the attribute is not an
// Epetra_Vector, throw a PythonException.
Teuchos::RCP<Epetra_Vector>
getEpetraVectorObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return a reference
// counted pointer to the const Epetra_Vector value of the attribute.
// If the attribute does not exist or the attribute is not an
// Epetra_Vector, throw a PythonException.
Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return a reference
// counted pointer to the const Epetra_Vector value of the i-th item
// of the attribute.  If the attribute does not exist or the attribute
// is not a sequence of Epetra_Vectors, throw a PythonException.
Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_MultiVector value of the attribute.
// If the attribute does not exist or the attribute is not an
// Epetra_MultiVector, throw a PythonException.
Teuchos::RCP<Epetra_MultiVector>
getEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return a reference
// counted pointer to the const Epetra_MultiVector value of the
// attribute.  If the attribute does not exist or the attribute is not
// an Epetra_MultiVector, throw a PythonException.
Teuchos::RCP<const Epetra_MultiVector>
getConstEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_Operator value of the attribute.  If
// the attribute does not exist or the attribute is not an
// Epetra_Operator, throw a PythonException.
Teuchos::RCP<Epetra_Operator>
getEpetraOperatorObjectAttr(PyObject * object, CONST char * name);

#endif  // HAVE_TEUCHOS

////////////////////////////////////////////////////////////////////////

}  // Namespace PyTrilinos

#endif // EPETRA_PYUTIL_H
