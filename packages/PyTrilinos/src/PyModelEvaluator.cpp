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

#include "PyModelEvaluator.h"
#include "swigpyrun.h"
#include <algorithm>

#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
// #include "Epetra_FastCrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
//#include "Epetra_MsrMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"

// PyModelEvaluator::PyModelEvaluator(PyObject* modelProperties) :
//   _modelProperties(NULL)
// {
//   // Obtain the class object for ModelProperties
//   PyObject * classModelProperties = getClassFromGlobals("ModelProperties");
//   if (!classModelProperties) throw PythonException();

//   // Ensure constructor argument is a ModelProperty object
//   if (!PyObject_IsInstance(modelProperties, classModelProperties))
//   {
//     PyErr_SetString(PyExc_TypeError,
// 		    "argument 'modelProperties' is not of type ModelProperties");
//     throw PythonException();
//   }

//   // Ensure the constructor argument has required attributes with non-None values
//   const char * required[ ] = {"x", "f", "evalModel", NULL};
//   int i = 0;
//   while(required[i] != NULL)
//   {
//     if (objectAttrIsNone(modelProperties, required[i]))
//     {
//       PyErr_Format(PyExc_ValueError, "attribute '%s' is required", required[i]);
//       throw PythonException();
//     }
//     ++i;
//   }

//   // Constructor argument is OK, so store it and increase its reference count
//   _modelProperties = modelProperties;
//   Py_INCREF(_modelProperties);
// }

// PyModelEvaluator::~PyModelEvaluator()
// {
//   // Release the stored python object
//   Py_XDECREF(_modelProperties);
// }

// Teuchos::RCP<const Epetra_Map>
// PyModelEvaluator::get_x_map() const
// {
//   Teuchos::RCP<const Epetra_Vector> x =
//     getConstEpetraVectorObjectAttr(_modelProperties, "x");
//   return getEpetraMapPtrFromEpetraBlockMap(x->Map());
// }

// Teuchos::RCP<const Epetra_Map>
// PyModelEvaluator::get_f_map() const
// {
//   EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f =
//     getEvaluationObjectAttr(_modelProperties, "f");
//   return getEpetraMapPtrFromEpetraBlockMap(f->Map());
// }

// Teuchos::RCP<const Epetra_Map>
// PyModelEvaluator::get_p_map(int l) const
// {
//   if (objectAttrIsNone(_modelProperties, "p"))
//     return EpetraExt::ModelEvaluator::get_p_map(l);
//   Teuchos::RCP<const Epetra_Vector> p =
//     getConstEpetraVectorItemObjectAttr(_modelProperties, "p", l);
//   return getEpetraMapPtrFromEpetraBlockMap(p->Map());
// }

// // Teuchos::RCP<const Teuchos::Array<std::string> >
// // PyModelEvaluator::get_p_names(int l) const
// // {
// // }

// Teuchos::RCP<const Epetra_Map>
// PyModelEvaluator::get_g_map(int j) const
// {
//   EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> g =
//     getEvaluationItemObjectAttr(_modelProperties, "g", j);
//   return getEpetraMapPtrFromEpetraBlockMap(g->Map());
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_x_init() const
// {
//   if (objectAttrIsNone(_modelProperties, "x_init"))
//     return EpetraExt::ModelEvaluator::get_x_init();
//   return getConstEpetraVectorObjectAttr(_modelProperties, "x_init");
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_x_dot_init() const
// {
//   if (objectAttrIsNone(_modelProperties, "x_dot_init"))
//     return EpetraExt::ModelEvaluator::get_x_dot_init();
//   return getConstEpetraVectorObjectAttr(_modelProperties, "x_dot_init");
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_p_init(int l) const
// {
//   if (objectAttrIsNone(_modelProperties, "p_init"))
//     return EpetraExt::ModelEvaluator::get_p_init(l);
//   return getConstEpetraVectorItemObjectAttr(_modelProperties, "p_init", l);
// }

// double
// PyModelEvaluator::get_t_init() const
// {
//   if (objectAttrIsNone(_modelProperties, "t_init"))
//     return EpetraExt::ModelEvaluator::get_t_init();
//   return getFloatObjectAttr(_modelProperties, "t_init");
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_x_lower_bounds() const
// {
//   if (objectAttrIsNone(_modelProperties, "x_lower_bounds"))
//     return EpetraExt::ModelEvaluator::get_x_lower_bounds();
//   return getConstEpetraVectorObjectAttr(_modelProperties, "x_lower_bounds");
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_x_upper_bounds() const
// {
//   if (objectAttrIsNone(_modelProperties, "x_upper_bounds"))
//     return EpetraExt::ModelEvaluator::get_x_upper_bounds();
//   return getConstEpetraVectorObjectAttr(_modelProperties, "x_upper_bounds");
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_p_lower_bounds(int l) const
// {
//   if (objectAttrIsNone(_modelProperties, "p_lower_bounds"))
//     return EpetraExt::ModelEvaluator::get_p_lower_bounds(l);
//   return getConstEpetraVectorItemObjectAttr(_modelProperties, "p_lower_bounds", l);
// }

// Teuchos::RCP<const Epetra_Vector>
// PyModelEvaluator::get_p_upper_bounds(int l) const
// {
//   if (objectAttrIsNone(_modelProperties, "p_upper_bounds"))
//     return EpetraExt::ModelEvaluator::get_p_upper_bounds(l);
//   return getConstEpetraVectorItemObjectAttr(_modelProperties, "p_upper_bounds", l);
// }

// double
// PyModelEvaluator::get_t_lower_bound() const
// {
//   if (objectAttrIsNone(_modelProperties, "t_lower_bound"))
//     return EpetraExt::ModelEvaluator::get_t_lower_bound();
//   return getFloatObjectAttr(_modelProperties, "t_lower_bound");
// }

// double
// PyModelEvaluator::get_t_upper_bound() const
// {
//   if (objectAttrIsNone(_modelProperties, "t_upper_bound"))
//     return EpetraExt::ModelEvaluator::get_t_upper_bound();
//   return getFloatObjectAttr(_modelProperties, "t_upper_bound");
// }

// Teuchos::RCP<Epetra_Operator>
// PyModelEvaluator::create_W() const
// {
//   if (objectAttrIsNone(_modelProperties, "W"))
//     return EpetraExt::ModelEvaluator::create_W();
//   return getEpetraOperatorObjectAttr(_modelProperties, "W");
// }

// Teuchos::RCP<Epetra_Operator>
// PyModelEvaluator::create_DfDp_op(int l) const
// {
//   if (objectAttrIsNone(_modelProperties, "DfDp"))
//     return EpetraExt::ModelEvaluator::create_DfDp_op(l);
//   return getDerivativeItemObjectAttr(_modelProperties, "DfDp", l).getLinearOp();
// }

// Teuchos::RCP<Epetra_Operator>
// PyModelEvaluator::create_DgDx_dot_op(int j) const
// {
//   if (objectAttrIsNone(_modelProperties, "DgDx_dot"))
//     return EpetraExt::ModelEvaluator::create_DgDx_dot_op(j);
//   return getDerivativeItemObjectAttr(_modelProperties, "DgDx_dot", j).getLinearOp();
// }

// Teuchos::RCP<Epetra_Operator>
// PyModelEvaluator::create_DgDx_op(int j) const
// {
//   if (objectAttrIsNone(_modelProperties, "DgDx"))
//     return EpetraExt::ModelEvaluator::create_DgDx_op(j);
//   return getDerivativeItemObjectAttr(_modelProperties, "DgDx", j).getLinearOp();
// }

// // Teuchos::RCP<Epetra_Operator>
// // PyModelEvaluator::create_DgDp_op( int j, int l ) const
// // {
// //   if (objectAttrIsNone(_modelProperties, "DgDp"))
// //     return EpetraExt::ModelEvaluator::create_DgDp_op();
// // }

// EpetraExt::ModelEvaluator::InArgs
// PyModelEvaluator::createInArgs() const
// {
//   InArgsSetup inArgs;
//   // x attribute
//   if (!objectAttrIsNone(_modelProperties, "x"))
//   {
//     inArgs.setSupports(IN_ARG_x, true);
//     inArgs.set_x(getConstEpetraVectorObjectAttr(_modelProperties, "x"));
//   }
//   // x_dot attribute
//   if (!objectAttrIsNone(_modelProperties, "x_dot"))
//   {
//     inArgs.setSupports(IN_ARG_x_dot, true);
//     inArgs.set_x(getConstEpetraVectorObjectAttr(_modelProperties, "x_dot"));
//   }
//   // p attribute
//   if (!objectAttrIsNone(_modelProperties, "p"))
//   {
//     PyObject * pObj = PyObject_GetAttrString(_modelProperties, "p");
//     if (!pObj) throw PythonException();
//     int Np = (int) PySequence_Length(pObj);
//     Py_DECREF(pObj);
//     if (Np < 0) throw PythonException();
//     inArgs.set_Np(Np);
//     for (int i=0; i < Np; ++i)
//       inArgs.set_p(i, getConstEpetraVectorItemObjectAttr(_modelProperties, "p", i));
//   }
//   // t attribute
//   if (!objectAttrIsNone(_modelProperties, "t"))
//   {
//     inArgs.setSupports(IN_ARG_t, true);
//     inArgs.set_t(getFloatObjectAttr(_modelProperties, "t"));
//   }
//   // alpha attribute
//   if (!objectAttrIsNone(_modelProperties, "alpha"))
//   {
//     inArgs.setSupports(IN_ARG_alpha, true);
//     inArgs.set_t(getFloatObjectAttr(_modelProperties, "alpha"));
//   }
//   // beta attribute
//   if (!objectAttrIsNone(_modelProperties, "beta"))
//   {
//     inArgs.setSupports(IN_ARG_beta, true);
//     inArgs.set_t(getFloatObjectAttr(_modelProperties, "beta"));
//   }
//   return inArgs;
// }

// EpetraExt::ModelEvaluator::OutArgs
// PyModelEvaluator::createOutArgs() const
// {
//   OutArgsSetup outArgs;
//   // Np and Ng attributes
//   int Np = 0;
//   int Ng = 0;
//   if (!objectAttrIsNone(_modelProperties, "p"))
//   {
//     PyObject * p = PyObject_GetAttrString(_modelProperties, "p");
//     if (!p) throw PythonException();
//     Np = (int) PySequence_Length(p);
//     Py_DECREF(p);
//     if (Np < 0) throw PythonException();
//   }
//   if (!objectAttrIsNone(_modelProperties, "g"))
//   {
//     PyObject * g = PyObject_GetAttrString(_modelProperties, "g");
//     if (!g) throw PythonException();
//     Ng = (int) PySequence_Length(g);
//     Py_DECREF(g);
//     if (Ng < 0) throw PythonException();
//   }
//   outArgs.set_Np_Ng(Np, Ng);
//   // f attribute
//   if (!objectAttrIsNone(_modelProperties, "f"))
//   {
//     outArgs.setSupports(OUT_ARG_f, true);
//     outArgs.set_f(getEvaluationObjectAttr(_modelProperties, "f"));
//   }
//   // W attributes
//   if (!objectAttrIsNone(_modelProperties, "W"))
//   {
//     outArgs.setSupports(OUT_ARG_W, true);
//     if (!objectAttrIsNone(_modelProperties, "W_properties"))
//       outArgs.set_W_properties(getDerivativePropertiesObjectAttr(_modelProperties,
// 								 "W_properties"));
//   }
//   // DfDp attributes
//   if (!objectAttrIsNone(_modelProperties, "DfDp_support"))
//   {
//     PyObject * DfDp_support = PyObject_GetAttrString(_modelProperties, "DfDp_support");
//     if (!DfDp_support) throw PythonException();
//     PyObject * DfDp_properties = PyObject_GetAttrString(_modelProperties, "DfDp_properties");
//     if (!DfDp_properties)
//     {
//       Py_DECREF(DfDp_support);
//       throw PythonException();
//     }
//     int n = (int) std::min(PySequence_Length(DfDp_support),
// 			   PySequence_Length(DfDp_properties));
//     Py_DECREF(DfDp_support);
//     Py_DECREF(DfDp_properties);
//     if (n < 0) throw PythonException();
//     for (int i=0; i < n; ++i)
//     {
//       outArgs.setSupports(OUT_ARG_DfDp, i,
// 			  getDerivativeSupportItemObjectAttr(_modelProperties,
// 							     "DfDp_support", i));
//       outArgs.set_DfDp_properties(i,
// 				  getDerivativePropertiesItemObjectAttr(_modelProperties,
// 									"DfDp_properties", i));
//     }
//   }
//   // DgDx attributes
//   if (!objectAttrIsNone(_modelProperties, "DgDx_support"))
//   {
//     PyObject * DgDx_support = PyObject_GetAttrString(_modelProperties, "DgDx_support");
//     if (!DgDx_support) throw PythonException();
//     PyObject * DgDx_properties = PyObject_GetAttrString(_modelProperties, "DgDx_properties");
//     if (!DgDx_properties)
//     {
//       Py_DECREF(DgDx_support);
//       throw PythonException();
//     }
//     int n = (int) std::min(PySequence_Length(DgDx_support),
// 			   PySequence_Length(DgDx_properties));
//     Py_DECREF(DgDx_support);
//     Py_DECREF(DgDx_properties);
//     if (n < 0) throw PythonException();
//     for (int i=0; i < n; ++i)
//     {
//       outArgs.setSupports(OUT_ARG_DgDx, i,
// 			  getDerivativeSupportItemObjectAttr(_modelProperties,
// 							     "DgDx_support", i));
//       outArgs.set_DgDx_properties(i,
// 				  getDerivativePropertiesItemObjectAttr(_modelProperties,
// 									"DgDx_properties", i));
//     }
//   }
//   // DgDx_dot attributes
//   if (!objectAttrIsNone(_modelProperties, "DgDx_dot_support"))
//   {
//     PyObject * DgDx_dot_support = PyObject_GetAttrString(_modelProperties,
// 							 "DgDx_dot_support");
//     if (!DgDx_dot_support) throw PythonException();
//     PyObject * DgDx_dot_properties = PyObject_GetAttrString(_modelProperties,
// 							    "DgDx_dot_properties");
//     if (!DgDx_dot_properties)
//     {
//       Py_DECREF(DgDx_dot_support);
//       throw PythonException();
//     }
//     int n = (int) std::min(PySequence_Length(DgDx_dot_support),
// 			   PySequence_Length(DgDx_dot_properties));
//     Py_DECREF(DgDx_dot_support);
//     Py_DECREF(DgDx_dot_properties);
//     if (n < 0) throw PythonException();
//     for (int i=0; i < n; ++i)
//     {
//       outArgs.setSupports(OUT_ARG_DgDx_dot, i,
// 			  getDerivativeSupportItemObjectAttr(_modelProperties,
// 							     "DgDx_dot_support", i));
//       outArgs.set_DgDx_dot_properties(
//         i, getDerivativePropertiesItemObjectAttr(_modelProperties,
// 						 "DgDx_dot_properties", i));
//     }
//   }
//   return outArgs;
// }

// void
// PyModelEvaluator::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
// {
//   // Convert inArgs and outArgs to a python argument list
//   PyObject * args = Py_BuildValue("(OO)",
// 				  convertInArgsToPython( inArgs ),
// 				  convertOutArgsToPython(outArgs));

//   // Extract and call the evalModel function
//   PyObject * evalModelObj = PyObject_GetAttrString(_modelProperties, "modelEval");
//   if (!evalModelObj)
//   {
//     Py_DECREF(args);
//     throw PythonException();
//   }
//   PyObject * result = PyObject_CallObject(evalModelObj, args);

//   // Exit: destruct evalModelObj and args; if result is NULL, then a
//   // python exception occurred calling modelEval(), so throw a
//   // PythonException; else result should be python None -- regardless,
//   // destruct it.
//   Py_DECREF(evalModelObj);
//   Py_DECREF(args);
//   if (!result) throw PythonException();
//   Py_DECREF(result);
// }

Teuchos::RCP<const Epetra_Map>
getEpetraMapPtrFromEpetraBlockMap(const Epetra_BlockMap & ebm)
{
  const Epetra_Map * em_ptr  = dynamic_cast<const Epetra_Map*>(&ebm);
  if (!em_ptr)
  {
    PyErr_SetString(PyExc_TypeError, "Cannot upcast BlockMap to Map");
    throw PythonException();
  }
  return Teuchos::rcp(em_ptr, false);
}

PyObject * getObjectFromGlobals(const char * name)
{
  PyObject * globals = PyEval_GetGlobals();
  if (!globals) return NULL;
  return PyDict_GetItemString(globals, name);
}

PyObject * getClassFromGlobals(const char * name)
{
  PyObject * object = getObjectFromGlobals(name);
  if (!object) return NULL;
  if (!PyType_Check(object))
  {
    PyErr_Format(PyExc_TypeError, "Object '%s' is not a class type", name);
    return NULL;
  }
  return object;
}

bool objectAttrIsNone(PyObject * object, const char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  bool result = (value == Py_None);
  Py_DECREF(value);
  return result;
}

bool objectAttrIsTrue(PyObject * object, const char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  bool result = (value == Py_True);
  Py_DECREF(value);
  return result;
}

bool getBoolObjectAttr(PyObject * object, const char * name)
{
  bool result;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyBool_Check(value))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type boolean", name);
    throw PythonException();
  }
  if (value == Py_True) result = true;
  else                  result = false;
  Py_DECREF(value);
  return result;
}

int getIntObjectAttr(PyObject * object, const char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  int result = (int) PyInt_AsLong(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

double getFloatObjectAttr(PyObject * object, const char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  double result = PyFloat_AsDouble(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

const char* getStringObjectAttr(PyObject * object, const char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  const char * result = PyString_AsString(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorObjectAttr(PyObject * object, const char * name)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Epetra_Vector *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_EV_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Vector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<const Epetra_Vector> result = 
    Teuchos::rcp(reinterpret_cast<Epetra_Vector *>(argp), false);
  Py_DECREF(value);
  return result;
}

Teuchos::RCP<Epetra_Vector>
getEpetraVectorObjectAttr(PyObject * object, const char * name)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Epetra_Vector *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_EV_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Vector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_Vector> result = 
    Teuchos::rcp(reinterpret_cast<Epetra_Vector *>(argp), false);
  Py_DECREF(value);
  return result;
}

Teuchos::RCP<Epetra_MultiVector>
getEpetraMultiVectorObjectAttr(PyObject * object, const char * name)
{
  static swig_type_info * swig_EMV_ptr = SWIG_TypeQuery("Epetra_MultiVector *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_EMV_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.MultiVector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_MultiVector> result =
    Teuchos::rcp(reinterpret_cast<Epetra_MultiVector *>(argp), false);
  Py_DECREF(value);
  return result;
}

Teuchos::RCP<const Epetra_MultiVector>
getConstEpetraMultiVectorObjectAttr(PyObject * object, const char * name)
{
  static swig_type_info * swig_EMV_ptr = SWIG_TypeQuery("Epetra_MultiVector *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_EMV_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.MultiVector", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<const Epetra_MultiVector> result =
    Teuchos::rcp(reinterpret_cast<Epetra_MultiVector *>(argp), false);
  Py_DECREF(value);
  return result;
}

Teuchos::RCP<Epetra_Operator>
getEpetraOperatorObjectAttr(PyObject * object, const char * name)
{
  static swig_type_info * swig_EO_ptr = SWIG_TypeQuery("Epetra_Operator *");
  void * argp;

  PyObject * value = PyObject_GetAttrString(object, name);
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_EO_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type Epetra.Operator", name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_Operator> result = 
    Teuchos::rcp(reinterpret_cast<Epetra_Operator *>(argp), false);
  Py_DECREF(value);
  return result;
}

EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationObjectAttr(PyObject * object, const char * name)
{
  // The Evaluation python object
  PyObject * classEvaluation = getClassFromGlobals("Evaluation");
  if (!classEvaluation) throw PythonException();
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
  const char * typeStr = getStringObjectAttr(value, "type");
  if (typeStr == "exact")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT;
  if (typeStr == "approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV;
  if (typeStr == "very_approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV;
  Py_DECREF(value);
  return EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>(vector, type);
}

EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportObjectAttr(PyObject * object, const char * name)
{
  PyObject * classDerivativeSupport = getClassFromGlobals("DerivativeSupport");
  if (!classDerivativeSupport) throw PythonException();
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

EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesObjectAttr(PyObject * object, const char * name)
{
  PyObject * classDerivativeProperties = getClassFromGlobals("DerivativeProperties");
  if (!classDerivativeProperties) throw PythonException();
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
  const char * linearity = getStringObjectAttr(value, "linearity");
  if (linearity == "unknown")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN;
  if (linearity == "const")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST;
  if (linearity == "nonconst")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST;
  // rank attribute
  const char * rank = getStringObjectAttr(value, "rank");
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

EpetraExt::ModelEvaluator::DerivativeMultiVector
getDerivativeMultiVectorObjectAttr(PyObject * object, const char * name)
{
  PyObject * classDerivativeMultiVector = getClassFromGlobals("DerivativeMultiVector");
  if (!classDerivativeMultiVector) throw PythonException();
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
  const char * linearity = getStringObjectAttr(value, "linearity");
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

EpetraExt::ModelEvaluator::Derivative
getDerivativeObjectAttr(PyObject * object, const char * name)
{
  PyObject * classDerivative = getClassFromGlobals("Derivative");
  if (!classDerivative) throw PythonException();
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

PyObject *
getTupleObjectAttr(PyObject * object, const char * name)
{
  PyObject * result = PyObject_GetAttrString(object, name);
  if (!PyTuple_Check(result))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type tuple", name);
    Py_DECREF(result);
    throw PythonException();
  }
  return result;
}

Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorItemObjectAttr(PyObject * object, const char * name, int i)
{
  static swig_type_info * swig_EV_ptr = SWIG_TypeQuery("Epetra_Vector *");
  void * argp;
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  if (!SWIG_CheckState(SWIG_Python_ConvertPtr(item, &argp, swig_EV_ptr, 0)))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not tuple of type Epetra.Vector", name);
    Py_DECREF(item);
    throw PythonException();
  }
  Py_DECREF(item);
  return Teuchos::rcp(reinterpret_cast<Epetra_Vector*>(argp), false);
}

const char *
getStringItemObjectAttr(PyObject * object, const char * name, int i)
{
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  const char * result = PyString_AsString(item);
  Py_DECREF(item);
  if (PyErr_Occurred()) throw PythonException();
  return result;
}

EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>
getEvaluationItemObjectAttr(PyObject * object, const char * name, int i)
{
  // The Evaluation python class object
  PyObject * classEvaluation = getClassFromGlobals("Evaluation");
  if (!classEvaluation) throw PythonException();
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
  const char * typeStr = getStringObjectAttr(item, "type");
  if (typeStr == "exact")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT;
  if (typeStr == "approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV;
  if (typeStr == "very_approx_deriv")
    type = EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV;
  Py_DECREF(item);
  return EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector>(vector, type);
}

EpetraExt::ModelEvaluator::DerivativeSupport
getDerivativeSupportItemObjectAttr(PyObject * object, const char * name, int i)
{
  // The DerivativeSupport python class object
  PyObject * classDerivativeSupport = getClassFromGlobals("DerivativeSupport");
  if (!classDerivativeSupport) throw PythonException();
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
  return result;}

EpetraExt::ModelEvaluator::DerivativeProperties
getDerivativePropertiesItemObjectAttr(PyObject * object, const char * name, int i)
{
  // The DerivativeProperties python class object
  PyObject * classDerivativeProperties = getClassFromGlobals("DerivativeProperties");
  if (!classDerivativeProperties) throw PythonException();
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
  const char * linearity = getStringObjectAttr(item, "linearity");
  if (linearity == "unknown")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN;
  if (linearity == "const")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST;
  if (linearity == "nonconst")
    result.linearity = EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST;
  // rank attribute
  const char * rank = getStringObjectAttr(item, "rank");
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

EpetraExt::ModelEvaluator::Derivative
getDerivativeItemObjectAttr(PyObject * object, const char * name, int i)
{
  // The Derivative python class object
  PyObject * classDerivative = getClassFromGlobals("Derivative");
  if (!classDerivative) throw PythonException();
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

PyObject *
convertEpetraMultiVectorToPython(const Epetra_MultiVector * emv)
{
  // SWIG initialization
  static swig_type_info * swig_EMV_ptr  = SWIG_TypeQuery("Epetra_MultiVector      *");
  static swig_type_info * swig_ENMV_ptr = SWIG_TypeQuery("Epetra_NumPyMultiVector *");

  const Epetra_NumPyMultiVector * enmv = dynamic_cast<const Epetra_NumPyMultiVector*>(emv);
  if (enmv) return SWIG_NewPointerObj((void*) enmv, swig_ENMV_ptr, 1);

  return SWIG_NewPointerObj((void*) emv, swig_EMV_ptr, 1);
}

PyObject *
convertEpetraVectorToPython(const Epetra_Vector * ev)
{
  // SWIG initialization
  static swig_type_info * swig_EV_ptr  = SWIG_TypeQuery("Epetra_Vector      *");
  static swig_type_info * swig_ENV_ptr = SWIG_TypeQuery("Epetra_NumPyVector *");

  const Epetra_NumPyVector * env = dynamic_cast<const Epetra_NumPyVector*>(ev);
  if (env) return SWIG_NewPointerObj((void*) env, swig_ENV_ptr, 1);

  return SWIG_NewPointerObj((void*) ev, swig_EV_ptr, 1);
}

PyObject *
convertEpetraOperatorToPython(Epetra_Operator * eo)
{
  // SWIG initialization
  static swig_type_info * swig_EO_ptr   = SWIG_TypeQuery("Epetra_Operator        *");
  //static swig_type_info * swig_EFCO_ptr = SWIG_TypeQuery("Epetra_FastCrsOperator *");
  static swig_type_info * swig_EIO_ptr  = SWIG_TypeQuery("Epetra_InvOperator     *");
  static swig_type_info * swig_ERM_ptr  = SWIG_TypeQuery("Epetra_RowMatrix       *");
  static swig_type_info * swig_EBRM_ptr = SWIG_TypeQuery("Epetra_BasicRowMatrix  *");
  static swig_type_info * swig_ECM_ptr  = SWIG_TypeQuery("Epetra_CrsMatrix       *");
  //static swig_type_info * swig_EMM_ptr  = SWIG_TypeQuery("Epetra_MsrMatrix       *");
  static swig_type_info * swig_EVM_ptr  = SWIG_TypeQuery("Epetra_VbrMatrix       *");
  static swig_type_info * swig_EVRM_ptr = SWIG_TypeQuery("Epetra_VbrRowMatrix    *");
  static swig_type_info * swig_EFVM_ptr = SWIG_TypeQuery("Epetra_FEVbrMatrix     *");
  static swig_type_info * swig_EFCM_ptr = SWIG_TypeQuery("Epetra_FECrsMatrix     *");
  static swig_type_info * swig_EJM_ptr  = SWIG_TypeQuery("Epetra_JadMatrix       *");

  Epetra_VbrRowMatrix * evrm = dynamic_cast<Epetra_VbrRowMatrix*>(eo);
  if (evrm) return SWIG_NewPointerObj((void*) evrm, swig_EVRM_ptr, 1);

  Epetra_FEVbrMatrix * efvm = dynamic_cast<Epetra_FEVbrMatrix*>(eo);
  if (efvm) return SWIG_NewPointerObj((void*) efvm, swig_EFVM_ptr, 1);

  Epetra_FECrsMatrix * efcm = dynamic_cast<Epetra_FECrsMatrix*>(eo);
  if (efcm) return SWIG_NewPointerObj((void*) efcm, swig_EFCM_ptr, 1);

  Epetra_JadMatrix * ejm = dynamic_cast<Epetra_JadMatrix*>(eo);
  if (ejm) return SWIG_NewPointerObj((void*) ejm, swig_EJM_ptr, 1);

  Epetra_BasicRowMatrix * ebrm = dynamic_cast<Epetra_BasicRowMatrix*>(eo);
  if (ebrm) return SWIG_NewPointerObj((void*) ebrm, swig_EBRM_ptr, 1);

  Epetra_CrsMatrix * ecm = dynamic_cast<Epetra_CrsMatrix*>(eo);
  if (ecm) return SWIG_NewPointerObj((void*) ecm, swig_ECM_ptr, 1);

  //Epetra_MsrMatrix * emm = dynamic_cast<Epetra_MsrMatrix*>(eo);
  //if (emm) return SWIG_NewPointerObj((void*) emm, swig_EMM_ptr, 1);

  Epetra_VbrMatrix * evm = dynamic_cast<Epetra_VbrMatrix*>(eo);
  if (evm) return SWIG_NewPointerObj((void*) evm, swig_EVM_ptr, 1);

  Epetra_RowMatrix * erm = dynamic_cast<Epetra_RowMatrix*>(eo);
  if (erm) return SWIG_NewPointerObj((void*) erm, swig_ERM_ptr, 1);

  Epetra_InvOperator * eio = dynamic_cast<Epetra_InvOperator*>(eo);
  if (eio) return SWIG_NewPointerObj((void*) eio, swig_EIO_ptr, 1);

  //Epetra_FastCrsOperator * efco = dynamic_cast<Epetra_FastCrsOperator*>(eo);
  //if (efco) return SWIG_NewPointerObj((void*) efco, swig_EFCO_ptr, 1);

  return SWIG_NewPointerObj((void*) eo, swig_EO_ptr, 1);
}

PyObject *
convertArrayOfIntToPython(const Teuchos::Array<int> & tai)
{
  PyObject * args = NULL;
  int        res  = 0;

  int size = tai.size();
  args = PyTuple_New(size);
  for (int i=0; i < size; ++i)
  {
    res = PyTuple_SetItem(args, Py_ssize_t(i), PyInt_FromLong(long(i)));
    if (res)
    {
      Py_DECREF(args);
      throw PythonException();
    }
  }
  PyObject * classTupleOfInt = getClassFromGlobals("tuple_of_int");
  if (!classTupleOfInt)
  {
    Py_DECREF(args);
    throw PythonException();
  }
  PyObject * result = PyInstance_New(classTupleOfInt,args,NULL);
  Py_DECREF(classTupleOfInt);
  Py_DECREF(args);
  if (!result) throw PythonException();

  return result;
}

PyObject *
convertInArgsToPython(const EpetraExt::ModelEvaluator::InArgs & inArgs)
{
  int res = 0;
  PyObject * obj;
  // Python class object
  PyObject * classInArgs = getClassFromGlobals("InArgs");
  if (!classInArgs) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * inArgsObj  = PyInstance_New(classInArgs, obj, NULL);
  Py_DECREF(classInArgs);
  Py_DECREF(obj);
  if (!inArgsObj) throw PythonException();

  // t attribute
  obj = PyFloat_FromDouble(inArgs.get_t());
  if (!obj) throw PythonException();
  res = PyObject_SetAttrString(inArgsObj, "t", obj);
  Py_DECREF(obj);
  if (res < 0) throw PythonException();

  // alpha attribute
  obj = PyFloat_FromDouble(inArgs.get_alpha());
  if (!obj) throw PythonException();
  res = PyObject_SetAttrString(inArgsObj, "alpha", obj);
  Py_DECREF(obj);
  if (res < 0) throw PythonException();

  // beta attribute
  obj = PyFloat_FromDouble(inArgs.get_beta());
  if (!obj) throw PythonException();
  res = PyObject_SetAttrString(inArgsObj, "beta", obj);
  Py_DECREF(obj);
  if (res < 0) throw PythonException();

  // x attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x))
  {
    Teuchos::RCP<const Epetra_Vector> xPtr = inArgs.get_x();
    res = PyObject_SetAttrString(inArgsObj, "x",
				 convertEpetraVectorToPython(xPtr.get()));
    if (res < 0) throw PythonException();
  }

  // x_dot attribute
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot))
  {
    Teuchos::RCP<const Epetra_Vector> x_dotPtr = inArgs.get_x_dot();
    res = PyObject_SetAttrString(inArgsObj, "x_dot",
				 convertEpetraVectorToPython(x_dotPtr.get()));
    if (res < 0) throw PythonException();
  }

  // p attribute
  int Np = inArgs.Np();
  if (Np > 0)
  {
    obj = PyTuple_New(Np);
    for (int i=0; i < Np; ++i)
    {
      Teuchos::RCP<const Epetra_Vector> pPtr = inArgs.get_p(i);
      res = PyTuple_SetItem(obj, Py_ssize_t(i), convertEpetraVectorToPython(pPtr.get()));
      if (res)
      {
	Py_DECREF(obj);
	throw PythonException();
      }
    }
    PyObject * classTupleOfVector = getClassFromGlobals("tuple_of_Vector");
    if (!classTupleOfVector)
    {
      Py_DECREF(obj);
      throw PythonException();
    }
    res = PyObject_SetAttrString(inArgsObj, "p", PyInstance_New(classTupleOfVector,obj,NULL));
    Py_DECREF(classTupleOfVector);
    Py_DECREF(obj);
    if (res) throw PythonException();
  }

  return inArgsObj;
}

PyObject *
convertEvaluationToPython(const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval)
{
  PyObject * obj;
  int res;
  // Python class object
  PyObject * classEvaluation = getClassFromGlobals("Evaluation");
  if (!classEvaluation) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * evalObj  = PyInstance_New(classEvaluation, obj, NULL);
  Py_DECREF(classEvaluation);
  Py_DECREF(obj);
  if (!evalObj) throw PythonException();

  // vector attribute
  res = PyObject_SetAttrString(evalObj, "vector", convertEpetraVectorToPython(eval.get()));
  if (res < 0)
  {
    Py_DECREF(evalObj);
    throw PythonException();
  }

  // type attribute
  obj = NULL;
  switch(eval.getType())
  {
  case EpetraExt::ModelEvaluator::EVAL_TYPE_EXACT:
    obj = PyString_FromString("exact"); break;
  case EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV:
    obj = PyString_FromString("approx_deriv"); break;
  case EpetraExt::ModelEvaluator::EVAL_TYPE_VERY_APPROX_DERIV:
    obj = PyString_FromString("very_approx_deriv");
  }
  res = PyObject_SetAttrString(evalObj, "type", obj);
  Py_XDECREF(obj);
  if (res < 0)
  {
    Py_DECREF(evalObj);
    throw PythonException();
  }

  return evalObj;
}

PyObject *
convertDerivativePropertiesToPython(const
				    EpetraExt::ModelEvaluator::DerivativeProperties & dProps)
{
  PyObject * obj = NULL;
  int        res = 0;
  // Python class object
  PyObject * classDerivativeProperties = getClassFromGlobals("DerivativeProperties");
  if (!classDerivativeProperties) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * dPropsObj  = PyInstance_New(classDerivativeProperties, obj, NULL);
  Py_DECREF(classDerivativeProperties);
  Py_DECREF(obj);
  if (!dPropsObj) throw PythonException();

  // linearity attribute
  obj = NULL;
  switch(dProps.linearity)
  {
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_UNKNOWN:
    obj = PyString_FromString("unknown"); break;
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_CONST:
    obj = PyString_FromString("const"); break;
  case EpetraExt::ModelEvaluator::DERIV_LINEARITY_NONCONST:
    obj = PyString_FromString("nonconst");
  }
  res = PyObject_SetAttrString(dPropsObj, "linearity", obj);
  Py_XDECREF(obj);
  if (res < 0)
  {
    Py_DECREF(dPropsObj);
    throw PythonException();
  }

  // rank attribute
  obj = NULL;
  switch(dProps.rank)
  {
  case EpetraExt::ModelEvaluator::DERIV_RANK_UNKNOWN:
    obj = PyString_FromString("unknown"); break;
  case EpetraExt::ModelEvaluator::DERIV_RANK_FULL:
    obj = PyString_FromString("full"); break;
  case EpetraExt::ModelEvaluator::DERIV_RANK_DEFICIENT:
    obj = PyString_FromString("deficient");
  }
  res = PyObject_SetAttrString(dPropsObj, "rank", obj);
  Py_XDECREF(obj);
  if (res < 0)
  {
    Py_DECREF(dPropsObj);
    throw PythonException();
  }

  // supportsAdjoint attribute
  obj = dProps.supportsAdjoint ? Py_True : Py_False;
  res = PyObject_SetAttrString(dPropsObj, "supportsAdjoint", obj);
  if (res < 0)
  {
    Py_DECREF(dPropsObj);
    throw PythonException();
  }

  return dPropsObj;
}

PyObject *
convertDerivativeMultiVectorToPython(const
				     EpetraExt::ModelEvaluator::DerivativeMultiVector & derivMV)
{
  int        res = 0;
  PyObject * obj = NULL;
  // Python class object
  PyObject * classDerivativeMultiVector = getClassFromGlobals("DerivativeMultiVector");
  if (!classDerivativeMultiVector) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * derivMVObj  = PyInstance_New(classDerivativeMultiVector, obj, NULL);
  Py_DECREF(classDerivativeMultiVector);
  Py_DECREF(obj);
  if (!derivMVObj) throw PythonException();

  // mutiVector attribute
  res = PyObject_SetAttrString(derivMVObj, "multiVector",
			       convertEpetraMultiVectorToPython(derivMV.getMultiVector().get()));
  if (res < 0)
  {
    Py_DECREF(derivMVObj);
    throw PythonException();
  }

  // orientation attribute
  switch(derivMV.getOrientation())
  {
  case EpetraExt::ModelEvaluator::DERIV_MV_BY_COL:
    obj = PyString_FromString("mv_by_col"); break;
  case EpetraExt::ModelEvaluator::DERIV_TRANS_MV_BY_ROW:
    obj = PyString_FromString("trans_mv_by_row");
  }
  res = PyObject_SetAttrString(derivMVObj, "orientation", obj);
  Py_XDECREF(obj);
  if (res < 0)
  {
    Py_DECREF(derivMVObj);
    throw PythonException();
  }

  // paramIndexes attribute
  res = PyObject_SetAttrString(derivMVObj, "paramIndexes",
			       convertArrayOfIntToPython(derivMV.getParamIndexes()));
  if (res < 0)
  {
    Py_DECREF(derivMVObj);
    throw PythonException();
  }

  return derivMVObj;
}

PyObject *
convertDerivativeToPython(const
			  EpetraExt::ModelEvaluator::Derivative & deriv)
{
  PyObject * obj = NULL;
  int        res = 0;
  // Python class object
  PyObject * classDerivative = getClassFromGlobals("Derivative");
  if (!classDerivative) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * derivObj  = PyInstance_New(classDerivative, obj, NULL);
  Py_DECREF(classDerivative);
  Py_DECREF(obj);
  if (!derivObj) throw PythonException();

  // operator attribute
  Epetra_Operator * eo = deriv.getLinearOp().get();
  if (eo)
  {
    res = PyObject_SetAttrString(derivObj, "operator", convertEpetraOperatorToPython(eo));
    if (res < 0)
    {
      Py_DECREF(derivObj);
      throw PythonException();
    }
  }

  // derivativeMultiVector attribute
  EpetraExt::ModelEvaluator::DerivativeMultiVector dmv = deriv.getDerivativeMultiVector();
  if (dmv.getMultiVector().get())
  {
    res = PyObject_SetAttrString(derivObj, "derivativeMultiVector",
				 convertDerivativeMultiVectorToPython(dmv));
    if (res < 0)
    {
      Py_DECREF(derivObj);
      throw PythonException();
    }
  }

  return derivObj;
}

PyObject *
convertOutArgsToPython(const EpetraExt::ModelEvaluator::OutArgs & outArgs)
{
  PyObject * obj = NULL;
  int        res = 0;
  // Python class object
  PyObject * classOutArgs = getClassFromGlobals("OutArgs");
  if (!classOutArgs) throw PythonException();

  // Create an instance
  obj = Py_BuildValue("()");
  PyObject * outArgsObj  = PyInstance_New(classOutArgs, obj, NULL);
  Py_DECREF(classOutArgs);
  Py_DECREF(obj);
  if (!outArgsObj) throw PythonException();

  // g attribute
  int Ng = outArgs.Ng();
  if (Ng > 0)
  {
    obj = PyTuple_New(Ng);
    for (int i=0; i < Ng; ++i)
    {
      const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval =
	outArgs.get_g(i);
      res = PyTuple_SetItem(obj, Py_ssize_t(i), convertEvaluationToPython(eval));
      if (res)
      {
	Py_DECREF(obj);
	throw PythonException();
      }
    }
    PyObject * classTupleOfEvaluation = getClassFromGlobals("tuple_of_Evaluation");
    if (!classTupleOfEvaluation)
    {
      Py_DECREF(obj);
      throw PythonException();
    }
    res = PyObject_SetAttrString(outArgsObj, "g",
				 PyInstance_New(classTupleOfEvaluation,obj,NULL));
    Py_DECREF(classTupleOfEvaluation);
    Py_DECREF(obj);
    if (res) throw PythonException();
  }

  // f attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_f))
  {
    const EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> eval =
      outArgs.get_f();
    res = PyObject_SetAttrString(outArgsObj, "f", convertEvaluationToPython(eval));
    if (res) throw PythonException();
  }

  // W attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W))
  {
    Teuchos::RCP<Epetra_Operator> eo = outArgs.get_W();
    res = PyObject_SetAttrString(outArgsObj, "W", convertEpetraOperatorToPython(eo.get()));
    if (res) throw PythonException();
  }

  // W_properties attribute
  if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W))
  {
    const EpetraExt::ModelEvaluator::DerivativeProperties W_properties =
      outArgs.get_W_properties();
    res = PyObject_SetAttrString(outArgsObj, "W_properties",
				 convertDerivativePropertiesToPython(W_properties));
    if (res) throw PythonException();
  }  

  return outArgsObj;
}

EpetraExt::ModelEvaluator::InArgs
convertInArgsFromPython(PyObject * source)
{
  EpetraExt::ModelEvaluator::InArgsSetup result;

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
    catch(PythonException &e)
    {
      PyErr_Clear();
    }
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
    catch(PythonException &e)
    {
      PyErr_Clear();
    }
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

EpetraExt::ModelEvaluator::OutArgs
convertOutArgsFromPython(PyObject * source)
{
  int Np = 0;
  int Ng = 0;
  EpetraExt::ModelEvaluator::OutArgsSetup result;

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
