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

#ifndef PYTRILINOS_EPETRAEXT_UTIL_HPP
#define PYTRILINOS_EPETRAEXT_UTIL_HPP

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

// Epetra includes
#include "Epetra_Vector.h"

// EpetraExt includes
#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_ModelEvaluator.h"

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////
// Helper functions.  Ultimately, we want to implement the functions
// 'convertInArgsFromPython', 'convertInArgsToPython',
// 'convertOutArgsFromPython' and 'convertOutArgsToPython' (where the
// '...ToPython' functions are friends of EpetraExt::ModelEvaluator).
// The following helper functions all support the imlpementation of
// these four conversion functions.
////////////////////////////////////////////////////////////////////////

// Define the string name of the PyTrilinos EpetraExt python module.
// This must be declared as non-const, because ultimately it gets
// passed to the PyImport_ImportModuleEx(...) function, which expects
// a non-const string (not that I know why...).
char PyTrilinosEpetraExt[21] = "PyTrilinos.EpetraExt";

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

}  //  Namespace PyTrilinos

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

namespace PyTrilinos
{

// Given a Teuchos Array of ints, return a tuple_of_int.  Return NULL
// upon failure.
PyObject *
convertArrayOfIntToPython(const Teuchos::Array<int> & tai);

}  // Namespace PyTrilinos

#endif // PYTRILINOS_EPETRAEXT_UTIL_HPP

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

