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
#ifdef HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"
#endif

// Epetra includes
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Operator.h"


////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

////////////////////////////////////////////////////////////

// Given an RCP to an Epetra_MultiVector, convert to a python object
// and return the pointer.  Attempt to downcast to an
// Epetra_NumPyMultiVector.
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< Epetra_MultiVector > *emv);

////////////////////////////////////////////////////////////

// Given an RCP to a const Epetra_MultiVector, convert to a python
// object and return the pointer.  Attempt to downcast to a const
// Epetra_NumPyMultiVector.
PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< const Epetra_MultiVector > *emv);

////////////////////////////////////////////////////////////

// Given an RCP to an Epetra_Vector, convert to a python object and
// return the pointer.  Attempt to downcast to an Epetra_NumPyVector.
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< Epetra_Vector > *ev);

////////////////////////////////////////////////////////////

// Given an RCP to a const Epetra_Vector, convert to a python object
// and return the pointer.  Attempt to downcast to a const
// Epetra_NumPyVector.
PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< const Epetra_Vector > *ev);

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

// Given a const Epetra_BlockMap &, return a reference counted pointer
// to a const Epetra_Map.  If the downcast cannot be performed, throw a
// PythonException.
Teuchos::RCP< const Epetra_Map >
getEpetraMapPtrFromEpetraBlockMap(const Epetra_BlockMap & ebm);

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
// counted pointer to the const Epetra_MultiVector value of the
// attribute.  If the attribute does not exist or the attribute is not
// an Epetra_MultiVector, throw a PythonException.
Teuchos::RCP< const Epetra_MultiVector >
getConstEpetraMultiVectorObjectAttr(PyObject   * object,
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

// Given a distributed array object, convert it to an RCP of a const
// Epetra_Map.
// Teuchos::RCP< const Epetra_Map >
// convertToEpetraMap(const Epetra_Comm & epetraComm,
//                    const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to an
// Epetra_MultiVector.
// Teuchos::RCP< Epetra_MultiVector >
// convertDistArrayToEpetraMultiVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to a
// const Epetra_MultiVector.
// Teuchos::RCP< const Epetra_MultiVector >
// convertDistArrayToConstEpetraMultiVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to an
// Epetra_Vector.
// Teuchos::RCP< Epetra_Vector >
// convertDistArrayToEpetraVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to a
// const Epetra_Vector.
// Teuchos::RCP< const Epetra_Vector >
// convertDistArrayToConstEpetraVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to an
// Epetra_FEVector.
// Teuchos::RCP< Epetra_FEVector >
// convertDistArrayToEpetraFEVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to a const
// Epetra_FEVector.
// Teuchos::RCP< const Epetra_FEVector >
// convertDistArrayToConstEpetraFEVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to an
// Epetra_IntVector.
// Teuchos::RCP< Epetra_IntVector >
// convertDistArrayToEpetraIntVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an RCP to a
// const Epetra_IntVector.
// Teuchos::RCP< const Epetra_IntVector >
// convertDistArrayToConstEpetraIntVector(PyObject * object);

////////////////////////////////////////////////////////////

#else   // HAVE_TEUCHOS

////////////////////////////////////////////////////////////

// Given a pointer to an Epetra_MultiVector, convert to a python
// object and return the pointer.  Attempt to downcast to an
// Epetra_NumPyMultiVector.
// PyObject *
// convertEpetraMultiVectorToDistArray(const Epetra_MultiVector * emv);

////////////////////////////////////////////////////////////

// Given a pointer to an Epetra_Vector, convert to a python object and
// return the pointer.  Attempt to downcast to an Epetra_NumPyVector.
// PyObject *
// convertEpetraVectorToDistArray(const Epetra_Vector * ev);

////////////////////////////////////////////////////////////

// Given a pointer to an Epetra_Operator, convert to a python object
// and return the pointer.  Attempt to downcast to any one of the many
// Epetra classes that derive from Epetra_Operator.
// PyObject *
// convertEpetraOperatorToDistArray(const Epetra_Operator * eo,
//                                  int cnvt_flags=0);

////////////////////////////////////////////////////////////

// Given a Python dimension data object, convert it to an Epetra_Map.
// const Epetra_Map *
// convertDimDataToEpetraMap(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an
// Epetra_MultiVector.
// Epetra_MultiVector *
// convertDistArrayToEpetraMultiVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an
// Epetra_Vector.
// Epetra_Vector *
// convertDistArrayToEpetraVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an
// Epetra_FEVector.
// Epetra_FEVector *
// convertDistArrayToEpetraFEVector(PyObject * object);

////////////////////////////////////////////////////////////

// Given a Python distributed array object, convert it to an
// Epetra_IntVector.
// Epetra_IntVector *
// convertDistArrayToEpetraIntVector(PyObject * object);

////////////////////////////////////////////////////////////

#endif  // HAVE_TEUCHOS

////////////////////////////////////////////////////////////

// Given a const Epetra_BlockMap &, return a Python dimension data
// object, which is a tuple of Python dimension data dictionaries that
// describe the Epetra_BlockMap, consistent with the DistArray
// Protocol.  The extraDim argument is to allow for the multiple
// vectors of an Epetra_MultiVector.  If an error occurs, return NULL.
// Note that an Epetra_BlockMap with variable element sizes is
// currently not supported and results in an error.
PyObject *
convertEpetraBlockMapToDimData(const Epetra_BlockMap & ebm,
                               int   extraDim=1);

////////////////////////////////////////////////////////////

// Given an Epetra_MultiVector, return a Python dictionary consistent
// with the DistArray Protocol.  If an error occurs, return NULL.
// PyObject *
// convertEpetraMultiVectorToDistArray(const Epetra_MultiVector & emv);

////////////////////////////////////////////////////////////////////////

}  // Namespace PyTrilinos

#endif // PYTRILINOS_EPETRA_UTIL_HPP
