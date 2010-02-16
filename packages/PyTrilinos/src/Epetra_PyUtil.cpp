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

// System includes
#include <algorithm>

// Epetra includes
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
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

// Local includes
#include "PyTrilinos_Util.h"
#include "Epetra_PyUtil.h"
#include "PythonException.h"
#include "swigpyrun.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"

////////////////////////////////////////////////////////////////////////

PyObject * convertEpetraMultiVectorToPython(const Epetra_MultiVector * emv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr = SWIG_TypeQuery("Epetra_NumPyMultiVector *");

  const Epetra_NumPyMultiVector * enmv = dynamic_cast<const Epetra_NumPyMultiVector*>(emv);
  if (enmv) enmv = new Epetra_NumPyMultiVector(View, *emv);
  return SWIG_NewPointerObj((void*) enmv, swig_ENMV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

PyObject * convertEpetraVectorToPython(const Epetra_Vector * ev)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr = SWIG_TypeQuery("Epetra_NumPyVector *");

  const Epetra_NumPyVector * env = dynamic_cast<const Epetra_NumPyVector*>(ev);
  if (!env) env = new Epetra_NumPyVector(View, *ev);
  return SWIG_NewPointerObj((void*) env, swig_ENV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

PyObject * convertEpetraOperatorToPython(Epetra_Operator * eo)
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

////////////////////////////////////////////////////////////////////////

#ifdef HAVE_TEUCHOS

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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Epetra_Vector>
getEpetraVectorObjectAttr(PyObject * object, CONST char * name)
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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorObjectAttr(PyObject * object, CONST char * name)
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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<const Epetra_Vector>
getConstEpetraVectorItemObjectAttr(PyObject * object, CONST char * name, int i)
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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Epetra_MultiVector>
getEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name)
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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<const Epetra_MultiVector>
getConstEpetraMultiVectorObjectAttr(PyObject * object, CONST char * name)
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

////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Epetra_Operator>
getEpetraOperatorObjectAttr(PyObject * object, CONST char * name)
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

////////////////////////////////////////////////////////////////////////

#endif   // HAVE_TEUCHOS
