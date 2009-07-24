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

#ifndef EPETRAEXT_PYUTIL_H
#define EPETRAEXT_PYUTIL_H

// Include the Python prototypes
#include "Python.h"

// Python developers improved the const-correctness of char* variables
// in the C-API with the advent of version 2.5.
#if PY_VERSION_HEX >= 0x02050000
#define CONST const
#else
#define CONST
#endif

// Teuchos includes
#include "Teuchos_RCPDecl.hpp"

// Epetra includes
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

// EpetraExt includes
#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_ModelEvaluator.h"

////////////////////////////////////////////////////////////////////////
// Helper functions.  Ultimately, we want to implement the two
// functions 'convertInArgsFromPython' and 'convertOutArgsFromPython'
// (which are friends of EpetraExt::ModelEvaluator).  The following
// helper functions all support the imlpementation of these two
// conversion functions.
////////////////////////////////////////////////////////////////////////

// Given the name of a python object, extract it from the global
// python namespace.  If it does not exist in the global namespace,
// return NULL.
PyObject * getObjectFromGlobals(CONST char * name);

// Given the name of a python class, extract it from the global
// python namespace.  If it does not exist in the global namespace,
// return NULL.
PyObject * getClassFromGlobals(CONST char * name);

// Given a python object and an attribute name, return 'true' if the
// attribute exists and is python 'None'.  If the attribute does not
// exist, throw a PythonException.
bool objectAttrIsNone(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return 'true' if the
// attribute exists and is python 'True'.  If the attribute does not
// exist, throw a PythonException.
bool objectAttrIsTrue(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the boolean
// value of the attribute.  If the attribute does not exist or the
// attribute is not a boolean, throw a PythonException.
bool getBoolObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the integer
// value of the attribute.  If the attribute does not exist or the
// attribute is not an integer, throw a PythonException.
int getIntObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the floating
// point (double) value of the attribute.  If the attribute does not
// exist or the attribute is not a float, throw a PythonException.
double getFloatObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the tuple value
// of the attribute.  If the attribute does not exist or the attribute
// is not a tuple, throw a PythonException.
PyObject * getTupleObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the string
// value of the attribute.  If the attribute does not exist or the
// attribute is not a string, throw a PythonException.
CONST char * getStringObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the string
// value of the i-th item of the attribute.  If the attribute does not
// exist or the attribute is not a sequence of strings, throw a
// PythonException.
CONST char * getStringItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a const Epetra_BlockMap, return a reference counted pointer
// to a const Epetra_Map.  If the upcast cannot be performed, throw a
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

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> value of the
// attribute.  If the attribute does not exist or the attribute is not
// an Evaluation, throw a PythonException.
EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> value of the
// i-th item of the attribute.  If the attribute does not exist or the
// attribute is not a sequence of Evaluations, throw a
// PythonException.
EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::DerivativeSupport value of the
// attribute.  If the attribute does not exist or the attribute is not
// a DerivativeSupport, throw a PythonException.
EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::DerivativeSupport value of the i-th item
// of the attribute.  If the attribute does not exist or the attribute
// is not a sequence of DerivativeSupports, throw a PythonException.
EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::DerivativeProperties value of the
// attribute.  If the attribute does not exist or the attribute is not
// a DerivativeProperties, throw a PythonException.
EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::DerivativeProperties value of the i-th
// item of the attribute.  If the attribute does not exist or the
// attribute is not a sequence of DerivativeProperties, throw a
// PythonException.
EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::DerivativeMultiVector value of the
// attribute.  If the attribute does not exist or the attribute is not
// a DerivativeMultiVector, throw a PythonException.
EpetraExt::ModelEvaluator::DerivativeMultiVector
getDerivativeMultiVectorObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::Derivative value of the attribute.  If
// the attribute does not exist or the attribute is not a Derivative,
// throw a PythonException.
EpetraExt::ModelEvaluator::Derivative 
getDerivativeObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the
// EpetraExt::ModelEvaluator::Derivative value of the attribute.  If
// the attribute does not exist or the attribute is not a sequence of
// Derivatives, throw a PythonException.
EpetraExt::ModelEvaluator::Derivative
getDerivativeItemObjectAttr(PyObject * object, CONST char * name, int i);

// Given a pointer to an Epetra_MultiVector, convert to a python
// object and return the pointer.  Attempt to downcast to an
// Epetra_NumPyMultiVector.
PyObject *
convertEpetraMultiVectorToPython(const Epetra_MultiVector * emv);

// Given a pointer to an Epetra_Vector, convert to a python object and
// return the pointer.  Attempt to downcast to an Epetra_NumPyVector.
PyObject * 
convertEpetraVectorToPython(const Epetra_Vector * ev);

// Given a pointer to an Epetra_Operator, convert to a python
// object and return the pointer.  Attempt to downcast to any one of
// the many Epetra classes that derive from Epetra_Operator.
PyObject *
convertEpetraOperatorToPython(Epetra_Operator * eo);

// Given a Teuchos Array of ints, return a tuple_of_int.  Return NULL
// upon failure.
PyObject *
convertArrayOfIntToPython(const Teuchos::Array<int> & tai);

// Given an EpetraExt::ModelEvaluator::InArgs object, convert to a
// python object and return the pointer.
PyObject *
convertInArgsToPython(const EpetraExt::ModelEvaluator::InArgs & inArgs);

// Given an EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
// object, convert to a python object and return the pointer.
PyObject *
convertEvaluationToPython(
    const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval);

// Given an EpetraExt::ModelEvaluator::DerivativeProperties object,
// convert to a python object and return the pointer.
PyObject *
convertDerivativePropertiesToPython(
    const EpetraExt::ModelEvaluator::DerivativeProperties & dProps);

// Given an EpetraExt::ModelEvaluator::DerivativeMultiVector object,
// convert to a python object and return the pointer.
PyObject *
convertDerivativeMultiVectorToPython(
    const EpetraExt::ModelEvaluator::DerivativeMultiVector & derivMV);

// Given an EpetraExt::ModelEvaluator::Derivative object, convert to a
// python object and return the pointer.
PyObject *
convertDerivativeToPython(
    const EpetraExt::ModelEvaluator::Derivative & deriv);

// Given an EpetraExt::ModelEvaluator::OutArgs object, convert to a
// python object and return the pointer.
PyObject *
convertOutArgsToPython(const EpetraExt::ModelEvaluator::OutArgs & outArgs);

// The following two conversion functions have already been declared
// as friends in the EpetraExt_ModelEvaluator.h header file.  However,
// some newer compilers require them to be declared explicitly again.
namespace EpetraExt
{
// Given a python InArgs object, return an equivalent
// EpetraExt::ModelEvaluator::InArgs object.  Upon failure, raise a
// PythonException.
ModelEvaluator::InArgs  convertInArgsFromPython( PyObject*);

// Given a python OutArgs object, return an equivalent
// EpetraExt::ModelEvaluator::OutArgs object.  Upon failure, raise a
// PythonException.
ModelEvaluator::OutArgs convertOutArgsFromPython(PyObject*);
}

#endif // EPETRAEXT_PYUTIL_H
