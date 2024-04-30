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

#ifndef PYTRILINOS_EPETRA_UTIL_HPP
#define PYTRILINOS_EPETRA_UTIL_HPP

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
#include "PyTrilinos_DAP.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Operator.h"

namespace PyTrilinos
{

////////////////////////////////////////////////////////////

// Given an RCP to an Epetra_MultiVector, convert to a python object
// and return the pointer.
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< Epetra_MultiVector > *emv);

////////////////////////////////////////////////////////////

// Given an RCP to a const Epetra_MultiVector, convert to a python
// object and return the pointer.
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< const Epetra_MultiVector > *emv);

////////////////////////////////////////////////////////////

// Given an RCP to an Epetra_Vector, convert to a python object and
// return the pointer.
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< Epetra_Vector > *ev);

////////////////////////////////////////////////////////////

// Given an RCP to a const Epetra_Vector, convert to a python object
// and return the pointer.
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< const Epetra_Vector > *ev);

////////////////////////////////////////////////////////////

// Attempt to convert a PyObject to an RCP to an Epetra_IntVector.
// The input PyObject could be a wrapped Epetra_IntVector, or a
// wrapped Domi::MDVector<int>, or an object that supports the
// DistArray Protocol, or, if the environment is serial, a simple
// NumPy array.
Teuchos::RCP< Epetra_IntVector > *
convertPythonToEpetraIntVector(PyObject * pyobj,
                               int * newmem);

////////////////////////////////////////////////////////////

// Attempt to convert a PyObject to an RCP to an Epetra_MultiVector.
// The input PyObject could be a wrapped Epetra_MultiVector, or a
// wrapped Domi::MDVector<double>, or an object that supports the
// DistArray Protocol, or, if the environment is serial, a simple
// NumPy array.
Teuchos::RCP< Epetra_MultiVector > *
convertPythonToEpetraMultiVector(PyObject * pyobj,
                                 int * newmem);

////////////////////////////////////////////////////////////

// Attempt to convert a PyObject to an RCP to an Epetra_Vector.  The
// input PyObject could be a wrapped Epetra_Vector, or a wrapped
// Domi::MDVector<double>, or an object that supports the DistArray
// Protocol, or, if the environment is serial, a simple NumPy array.
Teuchos::RCP< Epetra_Vector > *
convertPythonToEpetraVector(PyObject * pyobj,
                            int * newmem);

////////////////////////////////////////////////////////////

// Given an RCP to an Epetra_Operator, convert to a python object and
// return the pointer.  Attempt to downcast to any one of the many
// Epetra classes that derive from Epetra_Operator.
PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< Epetra_Operator > *eo);

////////////////////////////////////////////////////////////

// Given an RCP to a const Epetra_Operator, convert to a python object
// and return the pointer.  Attempt to downcast to any one of the many
// Epetra classes that derive from Epetra_Operator.
PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< const Epetra_Operator > *eo);

////////////////////////////////////////////////////////////

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_Vector value of the attribute.  If
// the attribute does not exist or the attribute is not an
// Epetra_Vector, throw a PythonException.
Teuchos::RCP< Epetra_Vector >
getEpetraVectorObjectAttr(PyObject   * object,
                          CONST char * name);

////////////////////////////////////////////////////////////

// Given a python object and an attribute name, return a reference
// counted pointer to the const Epetra_Vector value of the attribute.
// If the attribute does not exist or the attribute is not an
// Epetra_Vector, throw a PythonException.
Teuchos::RCP< const Epetra_Vector >
getConstEpetraVectorObjectAttr(PyObject   * object,
                               CONST char * name);

////////////////////////////////////////////////////////////

// Given a python object and an attribute name, return a reference
// counted pointer to the const Epetra_Vector value of the i-th item
// of the attribute.  If the attribute does not exist or the attribute
// is not a sequence of Epetra_Vectors, throw a PythonException.
Teuchos::RCP< const Epetra_Vector >
getConstEpetraVectorItemObjectAttr(PyObject   * object,
                                   CONST char * name,
                                   int          i);

////////////////////////////////////////////////////////////

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_MultiVector value of the attribute.
// If the attribute does not exist or the attribute is not an
// Epetra_MultiVector, throw a PythonException.
Teuchos::RCP< Epetra_MultiVector >
getEpetraMultiVectorObjectAttr(PyObject   * object,
                               CONST char * name);

////////////////////////////////////////////////////////////

// Given a python object and an attribute name, return a reference
// counted pointer to the Epetra_Operator value of the attribute.  If
// the attribute does not exist or the attribute is not an
// Epetra_Operator, throw a PythonException.
Teuchos::RCP< Epetra_Operator >
getEpetraOperatorObjectAttr(PyObject   * object,
                            CONST char * name);

////////////////////////////////////////////////////////////

// Given a const Epetra_BlockMap &, return a Python dimension data
// object, which is a tuple of Python dimension data dictionaries that
// describe the Epetra_BlockMap, consistent with the DistArray
// Protocol.  The extraDim argument is to allow for the multiple
// vectors of an Epetra_MultiVector.  If an error occurs, return NULL.
// Note that an Epetra_BlockMap with variable element sizes is
// currently not supported and results in an error.
PyObject *
convertToDimData(const Epetra_BlockMap & ebm,
                 int   extraDim=1);

////////////////////////////////////////////////////////////

// Given an Epetra_IntVector, return a Python dictionary consistent
// with the DistArray Protocol.  If an error occurs, return NULL.
PyObject *
convertToDistArray(const Epetra_IntVector & emv);

////////////////////////////////////////////////////////////

// Given an Epetra_MultiVector, return a Python dictionary consistent
// with the DistArray Protocol.  If an error occurs, return NULL.
PyObject *
convertToDistArray(const Epetra_MultiVector & emv);

////////////////////////////////////////////////////////////////////////

}  // Namespace PyTrilinos

#endif // PYTRILINOS_EPETRA_UTIL_HPP

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

