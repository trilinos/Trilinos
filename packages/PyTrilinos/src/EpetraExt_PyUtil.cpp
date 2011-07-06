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

// Local includes
#include "Epetra_PyUtil.h"
#include "PyTrilinos_config.h"
#include "PyTrilinos_Util.h"
#include "EpetraExt_PyUtil.h"
#include "PythonException.h"
#include "swigpyrun.h"

// Teuchos include
#ifdef HAVE_TEUCHOS
#include "Teuchos_Array.hpp"
#endif

// System includes
#include <algorithm>

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationObjectAttr(PyObject * object, CONST char * name)
{
  // The Evaluation python object
  static
  PyObject * classEvaluation = NULL;
  if (!classEvaluation)
  {
    classEvaluation = getClassFromModule(PyTrilinosEpetraExt, "Evaluation");
    if (!classEvaluation) throw PythonException();
  }
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyObject_IsInstance(value, classEvaluation))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Evaluation", name);
    Py_DECREF(value);
    throw PythonException();
  }
  // vector attribute
  Teuchos::RCP<Epetra_Vector> vector = getEpetraVectorObjectAttr(value, "vector");
  // type attribute
  EpetraExt::ModelEvaluator::EEvalType type;
  CONST char * typeStr = getStringObjectAttr(value, "type");
  if (typeStr == "exact")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT;
  if (typeStr == "approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV;
  if (typeStr == "very_approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV;
  Py_DECREF(value);
  return EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>(vector, type);
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  // The Evaluation python class object
  static
  PyObject * classEvaluation = NULL;
  if (!classEvaluation)
  {
    classEvaluation = getClassFromModule(PyTrilinosEpetraExt, "Evaluation");
    if (!classEvaluation) throw PythonException();
  }
  // Get the item from the object attribute
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  if (!PyObject_IsInstance(item, classEvaluation))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of Evaluation", name);
    Py_DECREF(item);
    throw PythonException();
  }
  // vector attribute
  Teuchos::RCP<Epetra_Vector> vector = getEpetraVectorObjectAttr(item, "vector");
  // type attribute
  EpetraExt::ModelEvaluator::EEvalType type;
  CONST char * typeStr = getStringObjectAttr(item, "type");
  if (typeStr == "exact")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT;
  if (typeStr == "approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV;
  if (typeStr == "very_approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV;
  Py_DECREF(item);
  return EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>(vector, type);
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportObjectAttr(PyObject * object, CONST char * name)
{
  static
  PyObject * classDerivativeSupport = NULL;
  if (!classDerivativeSupport)
  {
    classDerivativeSupport = getClassFromModule(PyTrilinosEpetraExt, "DerivativeSupport");
    if (!classDerivativeSupport) throw PythonException();
  }
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyObject_IsInstance(value, classDerivativeSupport))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type DerivativeSupport", name);
    Py_DECREF(value);
    throw PythonException();
  }
  EpetraExt::ModelEvaluator::EDerivativeLinearOp linearOp;
  EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation orientation;
  EpetraExt::ModelEvaluator::DerivativeSupport result;
  if (getBoolObjectAttr(value, "linearOp"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_LINEAR_OP);
  if (getBoolObjectAttr(value, "mVByCol"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_MV_BY_COL);
  if (getBoolObjectAttr(value, "transMVByRow"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  // The DerivativeSupport python class object
  static
  PyObject * classDerivativeSupport = NULL;
  if (!classDerivativeSupport)
  {
    classDerivativeSupport = getClassFromModule(PyTrilinosEpetraExt, "DerivativeSupport");
    if (!classDerivativeSupport) throw PythonException();
  }
  // Get the item from the object attribute
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  if (!PyObject_IsInstance(item, classDerivativeSupport))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of DerivativeSupport", name);
    Py_DECREF(item);
    throw PythonException();
  }
  EpetraExt::ModelEvaluator::EDerivativeLinearOp linearOp;
  EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation orientation;
  EpetraExt::ModelEvaluator::DerivativeSupport result;
  if (getBoolObjectAttr(item, "linearOp"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_LINEAR_OP);
  if (getBoolObjectAttr(item, "mVByCol"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_MV_BY_COL);
  if (getBoolObjectAttr(item, "transMVByRow"))
    result.plus(EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW);
  Py_DECREF(item);
  return result;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesObjectAttr(PyObject * object, CONST char * name)
{
  static
  PyObject * classDerivativeProperties = NULL;
  if (!classDerivativeProperties)
  {
    classDerivativeProperties = getClassFromModule(PyTrilinosEpetraExt, "DerivativeProperties");
    if (!classDerivativeProperties) throw PythonException();
  }
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyObject_IsInstance(value, classDerivativeProperties))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type DerivativeProperties", name);
    Py_DECREF(value);
    throw PythonException();
  }
  EpetraExt::ModelEvaluator::DerivativeProperties result;
  // linearity attribute
  CONST char * linearity = getStringObjectAttr(value, "linearity");
  if (linearity == "unknown")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN;
  if (linearity == "const")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST;
  if (linearity == "nonconst")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST;
  // rank attribute
  CONST char * rank = getStringObjectAttr(value, "rank");
  if (rank == "unknown")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN;
  if (rank == "full")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_FULL;
  if (rank == "deficient")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT;
  // supportsAdjoint attribute
  result.supportsAdjoint = getBoolObjectAttr(value, "supportsAdjoint");
  Py_DECREF(value);

  return result;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  // The DerivativeProperties python class object
  static
  PyObject * classDerivativeProperties = NULL;
  if (!classDerivativeProperties)
  {
    classDerivativeProperties = getClassFromModule(PyTrilinosEpetraExt, "DerivativeProperties");
    if (!classDerivativeProperties) throw PythonException();
  }
  // Get the item from the object attribute
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  if (!PyObject_IsInstance(item, classDerivativeProperties))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of DerivativeProperties", name);
    Py_DECREF(item);
    throw PythonException();
  }
  EpetraExt::ModelEvaluator::DerivativeProperties result;
  // linearity attribute
  CONST char * linearity = getStringObjectAttr(item, "linearity");
  if (linearity == "unknown")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN;
  if (linearity == "const")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST;
  if (linearity == "nonconst")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST;
  // rank attribute
  CONST char * rank = getStringObjectAttr(item, "rank");
  if (rank == "unknown")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN;
  if (rank == "full")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_FULL;
  if (rank == "deficient")
    result.rank = EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT;
  // supportsAdjoint attribute
  result.supportsAdjoint = getBoolObjectAttr(item, "supportsAdjoint");
  Py_DECREF(item);

  return result;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::DerivativeMultiVector
getDerivativeMultiVectorObjectAttr(PyObject * object, CONST char * name)
{
  static
  PyObject * classDerivativeMultiVector = NULL;
  if (!classDerivativeMultiVector)
  {
    classDerivativeMultiVector = getClassFromModule(PyTrilinosEpetraExt,
						    "DerivativeMultiVector");
    if (!classDerivativeMultiVector) throw PythonException();
  }
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyObject_IsInstance(value, classDerivativeMultiVector))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type DerivativeMultiVector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  // multiVector attribute
  Teuchos::RCP<Epetra_MultiVector> multiVector = 
    getEpetraMultiVectorObjectAttr(value, "multiVector");
  // orientation attribute
  EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation orientation;
  CONST char * linearity = getStringObjectAttr(value, "linearity");
  if (linearity == "mv_by_col")
    orientation = EpetraExt::ModelEvaluator::DERIV_MV_BY_COL;
  if (linearity == "trans_mv_by_row")
    orientation = EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW;
  // paramIndexes attribute
  PyObject * seq = PyObject_GetAttrString(value, "paramIndexes");
  if (!seq) throw PythonException();
  Py_ssize_t len = PySequence_Length(seq);
  if (len < 0) throw PythonException();
  Teuchos::Array<int> paramIndexes(len);
  for (Py_ssize_t i = 0; i < len; ++i)
  {
    PyObject * item = PySequence_GetItem(seq, i);
    if (!item) throw PythonException();
    paramIndexes[i] = (int) PyInt_AsLong(item);
    Py_DECREF(item);
    if (PyErr_Occurred()) throw PythonException();
  }
  Py_DECREF(seq);
  Py_DECREF(value);

  // Result
  return EpetraExt::ModelEvaluator::DerivativeMultiVector(multiVector, orientation,
							  paramIndexes);
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::Derivative
getDerivativeObjectAttr(PyObject * object, CONST char * name)
{
  static
  PyObject * classDerivative = NULL;
  if (!classDerivative)
  {
    classDerivative = getClassFromModule(PyTrilinosEpetraExt, "Derivative");
    if (!classDerivative) throw PythonException();
  }
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyObject_IsInstance(value, classDerivative))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Derivative", name);
    Py_DECREF(value);
    throw PythonException();
  }
  if (!objectAttrIsNone(value, "operator"))
  {
    // operator attribute
    Teuchos::RCP<Epetra_Operator> linearOp = 
      getEpetraOperatorObjectAttr(value, "operator");
    Py_DECREF(value);
    return EpetraExt::ModelEvaluator::Derivative(linearOp);
  }
  if (!objectAttrIsNone(value, "derivativeMultiVector"))
  {
    // derivativeMultiVector attribute
    EpetraExt::ModelEvaluator::DerivativeMultiVector derivativeMultiVector=
      getDerivativeMultiVectorObjectAttr(value, "derivativeMultiVector");
    Py_DECREF(value);
    return EpetraExt::ModelEvaluator::Derivative(derivativeMultiVector);
  }
  Py_DECREF(value);
  return EpetraExt::ModelEvaluator::Derivative();
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::Derivative
getDerivativeItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  // The Derivative python class object
  static
  PyObject * classDerivative = NULL;
  if (!classDerivative)
  {
    classDerivative = getClassFromModule(PyTrilinosEpetraExt, "Derivative");
    if (!classDerivative) throw PythonException();
  }
  // Get the item from the object attribute
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  if (!PyObject_IsInstance(item, classDerivative))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of Derivative", name);
    Py_DECREF(item);
    throw PythonException();
  }
  if (!objectAttrIsNone(item, "operator"))
  {
    // operator attribute
    Teuchos::RCP<Epetra_Operator> linearOp = 
      getEpetraOperatorObjectAttr(item, "operator");
    Py_DECREF(item);
    return EpetraExt::ModelEvaluator::Derivative(linearOp);
  }
  if (!objectAttrIsNone(item, "derivativeMultiVector"))
  {
    // derivativeMultiVector attribute
    EpetraExt::ModelEvaluator::DerivativeMultiVector derivativeMultiVector =
      getDerivativeMultiVectorObjectAttr(item, "derivativeMultiVector");
    Py_DECREF(item);
    return EpetraExt::ModelEvaluator::Derivative(derivativeMultiVector);
  }
  Py_DECREF(item);

  return EpetraExt::ModelEvaluator::Derivative();
}

////////////////////////////////////////////////////////////////////////

PyObject * convertInArgsToPython(const EpetraExt::ModelEvaluator::InArgs & inArgs)
{
  static
  PyObject * classInArgs 	= NULL;
  static
  PyObject * classTupleOfVector = NULL;
  PyObject * inArgsObj   	= NULL;
  PyObject * obj         	= NULL;
  PyObject * tupleOfVector      = NULL;
  int        res                = 0;
  int        Np                 = 0;
  Teuchos::RCP<const Epetra_Vector> xPtr;
  Teuchos::RCP<const Epetra_Vector> x_dotPtr;
  Teuchos::RCP<const Epetra_Vector> pPtr;

  // Python class objects
  if (!classInArgs)
  {
    classInArgs = getClassFromModule(PyTrilinosEpetraExt, "InArgs");
    if (!classInArgs) goto fail;
  }
  if (!classTupleOfVector)
  {
    classTupleOfVector = getClassFromModule(PyTrilinosEpetraExt, "tuple_of_Vector");
    if (!classTupleOfVector) goto fail;
  }

  // Create an instance
  inArgsObj = PyObject_CallObject(classInArgs, NULL);
  if (!inArgsObj) goto fail;

  // description attribute
  obj = PyString_FromString(inArgs.modelEvalDescription().c_str());
  if (!obj) goto fail;
  res = PyObject_SetAttrString(inArgsObj, "description", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  // t attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_t))
  {
    obj = PyFloat_FromDouble(inArgs.get_t());
    if (!obj) goto fail;
    res = PyObject_SetAttrString(inArgsObj, "t", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // alpha attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_alpha))
  {
    obj = PyFloat_FromDouble(inArgs.get_alpha());
    if (!obj) goto fail;
    res = PyObject_SetAttrString(inArgsObj, "alpha", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // beta attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_beta))
  {
    obj = PyFloat_FromDouble(inArgs.get_beta());
    if (!obj) goto fail;
    res = PyObject_SetAttrString(inArgsObj, "beta", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // x attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x))
  {
    xPtr = inArgs.get_x();
    obj  = convertEpetraVectorToPython(&xPtr);
    res  = PyObject_SetAttrString(inArgsObj, "x", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // x_dot attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot))
  {
    x_dotPtr = inArgs.get_x_dot();
    obj      = convertEpetraVectorToPython(&x_dotPtr);
    res      = PyObject_SetAttrString(inArgsObj, "x_dot", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // p attribute
  Np = inArgs.Np();
  if (Np > 0)
  {
    obj = PyTuple_New(Np);
    for (int i=0; i < Np; ++i)
    {
      pPtr = inArgs.get_p(i);
      res  = PyTuple_SetItem(obj, Py_ssize_t(i),
			     convertEpetraVectorToPython(&pPtr));
      if (res) goto fail;
    }
    tupleOfVector = PyObject_CallObject(classTupleOfVector, obj);
    res = PyObject_SetAttrString(inArgsObj, "p", tupleOfVector);
    Py_DECREF(tupleOfVector);
    Py_DECREF(obj);
    obj = NULL;
    if (res) goto fail;
  }

  return inArgsObj;

  fail:
  if (PyErr_Occurred()) PyErr_Print();
  Py_XDECREF(obj               );
  Py_XDECREF(inArgsObj         );
  Py_XDECREF(classTupleOfVector);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject * convertEvaluationToPython(
    const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval)
{
  static
  PyObject 	* classEvaluation = NULL;
  PyObject 	* obj             = NULL;
  PyObject 	* evalObj         = NULL;
  Epetra_Vector * vector          = NULL;
  int             res             = 0;
  Teuchos::RCP< Epetra_Vector > *smartarg = NULL;

  // Python class object
  if (!classEvaluation)
  {
    classEvaluation = getClassFromModule(PyTrilinosEpetraExt, "Evaluation");
    if (!classEvaluation) goto fail;
  }

  // Create an instance
  evalObj = PyObject_CallObject(classEvaluation, NULL);
  if (!evalObj) goto fail;

  // vector attribute
  vector = eval.get();
  if (vector)
  {
    smartarg = new Teuchos::RCP< Epetra_Vector>(vector, false);
    obj = convertEpetraVectorToPython(smartarg);
    res = PyObject_SetAttrString(evalObj, "vector", obj);
    Py_DECREF(obj);
    obj = NULL;
    delete smartarg;
    if (res < 0) goto fail;
  }

  // type attribute
  switch(eval.getType())
  {
  case EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT:
    obj = PyString_FromString("exact"); break;
  case EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV:
    obj = PyString_FromString("approx_deriv"); break;
  case EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV:
    obj = PyString_FromString("very_approx_deriv");
  }
  if (!obj) goto fail;
  res = PyObject_SetAttrString(evalObj, "type", obj);
  Py_XDECREF(obj);
  if (res < 0) goto fail;

  return evalObj;

  fail:
  Py_XDECREF(obj);
  Py_XDECREF(evalObj);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject * convertDerivativePropertiesToPython(
    const EpetraExt::ModelEvaluator::DerivativeProperties & dProps)
{
  static
  PyObject * classDerivativeProperties = NULL;
  PyObject * obj                       = NULL;
  PyObject * dPropsObj                 = NULL;
  int        res                       = 0;

  // Python class object
  if (!classDerivativeProperties)
  {
    classDerivativeProperties = getClassFromModule(PyTrilinosEpetraExt,
						   "DerivativeProperties");
    if (!classDerivativeProperties) goto fail;
  }

  // Create an instance
  dPropsObj = PyObject_CallObject(classDerivativeProperties, NULL);
  if (!dPropsObj) goto fail;

  // linearity attribute
  switch(dProps.linearity)
  {
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN:
    obj = PyString_FromString("unknown"); break;
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST:
    obj = PyString_FromString("const"); break;
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST:
    obj = PyString_FromString("nonconst");
  }
  if (!obj) goto fail;
  res = PyObject_SetAttrString(dPropsObj, "linearity", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  // rank attribute
  switch(dProps.rank)
  {
  case EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN:
    obj = PyString_FromString("unknown"); break;
  case EpetraExt::ModelEvaluator::DERIV_RANK_FULL:
    obj = PyString_FromString("full"); break;
  case EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT:
    obj = PyString_FromString("deficient");
  }
  if (!obj) goto fail;
  res = PyObject_SetAttrString(dPropsObj, "rank", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  // supportsAdjoint attribute
  obj = dProps.supportsAdjoint ? Py_True : Py_False;
  res = PyObject_SetAttrString(dPropsObj, "supportsAdjoint", obj);
  obj = NULL;
  if (res < 0) goto fail;

  return dPropsObj;

  fail:
  Py_XDECREF(obj      );
  Py_XDECREF(dPropsObj);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject * convertDerivativeMultiVectorToPython(
    const EpetraExt::ModelEvaluator::DerivativeMultiVector & derivMV)
{
  static
  PyObject * classDerivativeMultiVector = NULL;
  PyObject * obj                        = NULL;
  PyObject * derivMVObj                 = NULL;
  int        res                        = 0;
  Teuchos::RCP< Epetra_MultiVector > *smartarg = NULL;

  // Python class object
  if (!classDerivativeMultiVector)
  {
    classDerivativeMultiVector = getClassFromModule(PyTrilinosEpetraExt,
						    "DerivativeMultiVector");
    if (!classDerivativeMultiVector) goto fail;
  }

  // Create an instance
  derivMVObj = PyObject_CallObject(classDerivativeMultiVector, NULL);
  if (!derivMVObj) goto fail;

  // mutiVector attribute
  smartarg = new Teuchos::RCP< Epetra_MultiVector >(derivMV.getMultiVector());
  obj = convertEpetraMultiVectorToPython(smartarg);
  res = PyObject_SetAttrString(derivMVObj, "multiVector", obj);
  Py_DECREF(obj);
  obj = NULL;
  delete smartarg;
  if (res < 0) goto fail;

  // orientation attribute
  switch(derivMV.getOrientation())
  {
  case EpetraExt::ModelEvaluator::DERIV_MV_BY_COL:
    obj = PyString_FromString("mv_by_col"); break;
  case EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW:
    obj = PyString_FromString("trans_mv_by_row");
  }
  if (!obj) goto fail;
  res = PyObject_SetAttrString(derivMVObj, "orientation", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  // paramIndexes attribute
  obj = convertArrayOfIntToPython(derivMV.getParamIndexes());
  res = PyObject_SetAttrString(derivMVObj, "paramIndexes", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  return derivMVObj;

  fail:
  Py_XDECREF(obj       );
  Py_XDECREF(derivMVObj);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertDerivativeToPython(const EpetraExt::ModelEvaluator::Derivative & deriv)
{
  static
  PyObject * classDerivative = NULL;
  PyObject * obj             = NULL;
  PyObject * derivObj        = NULL;
  int        res             = 0;
  Teuchos::RCP< Epetra_Operator > eo;
  EpetraExt::ModelEvaluator::DerivativeMultiVector dmv;

  // Python class object
  if (!classDerivative)
  {
    classDerivative = getClassFromModule(PyTrilinosEpetraExt, "Derivative");
    if (!classDerivative) goto fail;
  }
  // Create an instance
  derivObj = PyObject_CallObject(classDerivative, NULL);
  if (!derivObj) goto fail;

  // operator attribute
  eo = deriv.getLinearOp();
  if (!eo.is_null())
  {
    obj = convertEpetraOperatorToPython(&eo);
    res = PyObject_SetAttrString(derivObj, "operator", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  // derivativeMultiVector attribute
  dmv = deriv.getDerivativeMultiVector();
  if (dmv.getMultiVector().get())
  {
    obj = convertDerivativeMultiVectorToPython(dmv);
    res = PyObject_SetAttrString(derivObj, "derivativeMultiVector", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res < 0) goto fail;
  }

  return derivObj;

  fail:
  Py_XDECREF(obj);
  Py_XDECREF(derivObj);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject * convertOutArgsToPython(const EpetraExt::ModelEvaluator::OutArgs & outArgs)
{
  static
  PyObject * classOutArgs 	    = NULL;
  static
  PyObject * classTupleOfEvaluation = NULL;
  PyObject * obj          	    = NULL;
  PyObject * tupleOfEval            = NULL;
  PyObject * outArgsObj   	    = NULL;
  int        res	            = 0;
  int        Ng 	            = 0;
  Teuchos::RCP< Epetra_Operator > WPtr;

  // Python class objects
  if (!classOutArgs)
  {
    classOutArgs = getClassFromModule(PyTrilinosEpetraExt, "OutArgs");
    if (!classOutArgs) goto fail;
  }
  if (!classTupleOfEvaluation)
  {
    classTupleOfEvaluation = getClassFromModule(PyTrilinosEpetraExt,
						"tuple_of_Evaluation");
    if (!classTupleOfEvaluation) goto fail;
  }

  // Create an instance
  outArgsObj = PyObject_CallObject(classOutArgs, NULL);
  if (!outArgsObj) goto fail;

  // description attribute
  obj = PyString_FromString(outArgs.modelEvalDescription().c_str());
  if (!obj) goto fail;
  res = PyObject_SetAttrString(outArgsObj, "description", obj);
  Py_DECREF(obj);
  obj = NULL;
  if (res < 0) goto fail;

  // g attribute
  Ng = outArgs.Ng();
  if (Ng > 0)
  {
    obj = PyTuple_New(Ng);
    for (int i=0; i < Ng; ++i)
    {
      res = PyTuple_SetItem(obj, Py_ssize_t(i),
			    convertEvaluationToPython(outArgs.get_g(i)));
      if (res) goto fail;
    }
    tupleOfEval = PyObject_CallObject(classTupleOfEvaluation, obj);
    res = PyObject_SetAttrString(outArgsObj, "g", tupleOfEval);
    Py_DECREF(tupleOfEval);
    Py_DECREF(obj);
    obj = NULL;
    if (res) goto fail;
  }

  // f attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_f))
  {
    obj = convertEvaluationToPython(outArgs.get_f());
    res = PyObject_SetAttrString(outArgsObj, "f", obj);
    Py_DECREF(obj);
    obj = NULL;
    if (res) goto fail;
  }

  // W attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W))
  {
    WPtr = outArgs.get_W();
    if (!WPtr.is_null())
    {
      obj = convertEpetraOperatorToPython(&WPtr);
      res = PyObject_SetAttrString(outArgsObj, "W", obj);
      Py_DECREF(obj);
      obj = NULL;
      if (res) goto fail;
    }
  }

  // W_properties attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W))
  {
    if (!WPtr.is_null())
    {
      obj = convertDerivativePropertiesToPython(outArgs.get_W_properties());
      res = PyObject_SetAttrString(outArgsObj, "W_properties", obj);
      Py_DECREF(obj);
      obj = NULL;
      if (res) goto fail;
    }
  }  

  return outArgsObj;

  fail:
  Py_XDECREF(obj       );
  Py_XDECREF(outArgsObj);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::InArgs
EpetraExt::convertInArgsFromPython(PyObject * source)
{
  EpetraExt::ModelEvaluator::InArgsSetup result;

  // description attribute
  result.setModelEvalDescription(std::string(getStringObjectAttr(source, "description")));

  // x attribute
  if (objectAttrIsTrue(source, "x"))
    result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_x, true);
  else
  {
    try
    {
      result.set_x(getConstEpetraVectorObjectAttr(source, "x"));
      result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_x, true);
    }
    catch(PythonException &e) { }
  }

  // x_dot attribute
  if (objectAttrIsTrue(source, "x_dot"))
    result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_x_dot, true);
  else
  {
    try
    {
      result.set_x(getConstEpetraVectorObjectAttr(source, "x_dot"));
      result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_x_dot, true);
    }
    catch(PythonException &e) { }
  }

  // p attribute
  PyObject * pObj = PyObject_GetAttrString(source, "p");
  if (!pObj) throw PythonException();
  if (PyInt_Check(pObj))
  {
    result.set_Np((int)PyInt_AsLong(pObj));
    Py_DECREF(pObj);
  }
  else
  {
    int Np = (int) PySequence_Length(pObj);
    Py_DECREF(pObj);
    if (Np < 0) throw PythonException();
    result.set_Np(Np);
    for (int i=0; i < Np; ++i)
      result.set_p(i, getConstEpetraVectorItemObjectAttr(source, "p", i));
  }

  // t attribute
  if (!objectAttrIsNone(source, "t"))
  {
    result.set_t(getFloatObjectAttr(source, "t"));
    result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_t, true);
  }

  // alpha attribute
  if (!objectAttrIsNone(source, "alpha"))
  {
    result.set_t(getFloatObjectAttr(source, "alpha"));
    result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_alpha, true);
  }

  // beta attribute
  if (!objectAttrIsNone(source, "beta"))
  {
    result.set_t(getFloatObjectAttr(source, "beta"));
    result.setSupports(EpetraExt::ModelEvaluator::IN_ARG_beta, true);
  }

  return result;
}

////////////////////////////////////////////////////////////////////////

EpetraExt::ModelEvaluator::OutArgs
EpetraExt::convertOutArgsFromPython(PyObject * source)
{
  int Np = 0;
  int Ng = 0;
  EpetraExt::ModelEvaluator::OutArgsSetup result;

  // description attribute
  result.setModelEvalDescription(std::string(getStringObjectAttr(source, "description")));

  // Number of p: Np
  PyObject * DfDpObj = PyObject_GetAttrString(source, "DfDp");
  if (!DfDpObj) throw PythonException();
  if (PyInt_Check(DfDpObj))
  {
    Np = (int) PyInt_AsLong(DfDpObj);
    Py_DECREF(DfDpObj);
  }
  else
  {
    Np = (int) PySequence_Length(DfDpObj);
    Py_DECREF(DfDpObj);
    if (Np < 0) throw PythonException();
  }

  // Number of g: Ng
  PyObject * gObj = PyObject_GetAttrString(source, "g");
  if (!gObj) throw PythonException();
  if (PyInt_Check(gObj))
  {
    Ng = (int) PyInt_AsLong(gObj);
    Py_DECREF(gObj);
  }
  else
  {
    Ng = (int) PySequence_Length(gObj);
    Py_DECREF(gObj);
    if (Ng < 0) throw PythonException();
  }

  // Size attributes
  result.set_Np_Ng(Np, Ng);

  // g attribute
  gObj = PyObject_GetAttrString(source, "g");
  if (!gObj) throw PythonException();
  bool dataProvided = (!PyInt_Check(gObj));
  Py_DECREF(gObj);
  if (dataProvided)
  {
    for (int i=0; i < Ng; ++i)
      result.set_g(i, getEvaluationItemObjectAttr(source, "g", i));
  }

  // f attribute
  if (objectAttrIsTrue(source, "f"))
    result.setSupports(EpetraExt::ModelEvaluator::OUT_ARG_f, true);
  else
  {
    try
    {
      result.set_f(getEvaluationObjectAttr(source, "f"));
      result.setSupports(EpetraExt::ModelEvaluator::OUT_ARG_f, true);
    }
    catch(PythonException &e)
    {
      PyErr_Clear();
    }
  }

  // W attribute
  if (objectAttrIsTrue(source, "W"))
    result.setSupports(EpetraExt::ModelEvaluator::OUT_ARG_W, true);
  else
  {
    try
    {
      result.set_f(getEvaluationObjectAttr(source, "W"));
      result.setSupports(EpetraExt::ModelEvaluator::OUT_ARG_W, true);
    }
    catch(PythonException &e)
    {
      PyErr_Clear();
    }
  }

  // DfDp attribute
  DfDpObj = PyObject_GetAttrString(source, "DfDp");
  if (!DfDpObj) throw PythonException();
  dataProvided = (!PyInt_Check(DfDpObj));
  Py_DECREF(DfDpObj);
  if (dataProvided)
  {
    for (int i=0; i < Np; ++i)
    {
      result.set_DfDp(i, getDerivativeItemObjectAttr(source, "DfDp", i));
      result.set_DfDp_properties(i, getDerivativePropertiesItemObjectAttr(source,
									  "DfDp_properties", i));
    }
  }

  // DgDx attribute
  PyObject * DgDxObj = PyObject_GetAttrString(source, "DgDx");
  if (!DgDxObj) throw PythonException();
  dataProvided = (!PyInt_Check(DgDxObj));
  Py_DECREF(DgDxObj);
  if (dataProvided)
  {
    for (int i=0; i < Ng; ++i)
    {
      result.set_DgDx(i, getDerivativeItemObjectAttr(source, "DgDx", i));
      result.set_DgDx_properties(i, getDerivativePropertiesItemObjectAttr(source,
							   "DgDx_properties", i));
    }
  }

  // DgDx_dot attribute
  PyObject * DgDx_dotObj = PyObject_GetAttrString(source, "DgDx_dot");
  if (!DgDx_dotObj) throw PythonException();
  dataProvided = (!PyInt_Check(DgDx_dotObj));
  Py_DECREF(DgDx_dotObj);
  if (dataProvided)
  {
    for (int i=0; i < Ng; ++i)
    {
      result.set_DgDx_dot(i, getDerivativeItemObjectAttr(source, "DgDx_dot", i));
      result.set_DgDx_dot_properties(i, getDerivativePropertiesItemObjectAttr(source,
							   "DgDx_dot_properties", i));
    }
  }

  // DgDp attribute
  // PyObject * DgDpObj = PyObject_GetAttrString(source, "DgDp");
  // if (!DgDpObj) throw PythonException();
  // bool dataProvided = (!PyInt_Check(DgDpObj));
  // Py_DECREF(DgDpObj);
  // if (dataProvided)
  // {
  //   for (int i=0; i < Ng; ++i)
  //   {
  //     result.set_DgDp(i, getDerivativeItemObjectAttr(source, "DgDp", i));
  //     result.set_DgDp_properties(i, getDerivativePropertiesItemObjectAttr(source,
  //							   "DgDp_properties", i));
  //   }
  // }

  return result;
}

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

PyObject * convertArrayOfIntToPython(const Teuchos::Array<int> & tai)
{
  static
  PyObject * classTupleOfInt = NULL;
  PyObject * args            = NULL;
  PyObject * result          = NULL;
  int        res             = 0;
  int        size            = 0;

  if (!classTupleOfInt)
  {
    classTupleOfInt = getClassFromModule(PyTrilinosEpetraExt, "tuple_of_int");
    if (!classTupleOfInt) goto fail;
  }
  size = tai.size();
  args = PyTuple_New(size);
  for (int i=0; i < size; ++i)
  {
    res = PyTuple_SetItem(args, Py_ssize_t(i), PyInt_FromLong(long(tai[i])));
    if (res) goto fail;
  }
  result = PyObject_CallObject(classTupleOfInt, args);
  Py_DECREF(args);
  args = NULL;
  if (!result) goto fail;

  return result;

  fail:
  Py_XDECREF(args);
  return NULL;
}

#endif   // HAVE_TEUCHOS
