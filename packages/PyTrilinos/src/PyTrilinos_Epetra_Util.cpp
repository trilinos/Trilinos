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

// Put local includes first to make sure that Python.h is the first included
// system header. This avoids compiler warnings about redefined variables
// such as _POSIX_C_SOURCE.
// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include "PyTrilinos_config.h"
#include "PyTrilinos_Util.hpp"
#include "PyTrilinos_Epetra_Util.hpp"
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_DAP.hpp"
#include "swigpyrun.h"

// System includes
#include <algorithm>

// Teuchos includes
#include "Teuchos_DefaultComm.hpp"

// Epetra includes
#include "Epetra_SerialComm.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"

#ifdef HAVE_PYTRILINOS_DOMI
#include "Domi_MDVector.hpp"
#include "PyTrilinos_Domi_Util.hpp"
#endif

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< Epetra_MultiVector > *emv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector >*");
  return SWIG_NewPointerObj((void*)emv, swig_ENMV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraMultiVectorToPython(const Teuchos::RCP< const Epetra_MultiVector > *cemv)
{
  // SWIG initialization
  static swig_type_info * swig_ENMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector >*");
  return SWIG_NewPointerObj((void*)cemv, swig_ENMV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< Epetra_Vector > *ev)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector >*");
  //
  // Convert to PyTrilinos::Epetra_Vector
  const Teuchos::RCP< Epetra_Vector > *env = new
    Teuchos::RCP< Epetra_Vector >(new Epetra_Vector(View, **ev, 0));
  return SWIG_NewPointerObj((void*)env, swig_ENV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraVectorToPython(const Teuchos::RCP< const Epetra_Vector > *cev)
{
  // SWIG initialization
  static swig_type_info * swig_ENV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector >*");
  //
  // Convert to PyTrilinos::Epetra_Vector
  const Teuchos::RCP< Epetra_Vector > *env = new
    Teuchos::RCP< Epetra_Vector >(new Epetra_Vector(View, **cev, 0));
  return SWIG_NewPointerObj((void*)env, swig_ENV_ptr, 1);
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_IntVector > *
convertPythonToEpetraIntVector(PyObject * pyobj,
                               int * newmem)
{
  // SWIG initialization
  static swig_type_info * swig_EIV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_IntVector >*");
  static swig_type_info * swig_DMDV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Domi::MDVector<int> >*");
  //
  // Result objects
  void *argp = 0;
  PyObject * distarray = 0;
  Teuchos::RCP< Epetra_IntVector > smartresult;
  Teuchos::RCP< Epetra_IntVector > * result;
#ifdef HAVE_PYTRILINOS_DOMI
  Teuchos::RCP< Domi::MDVector<int> > dmdv_rcp;
#endif
  *newmem = 0;
  //
  // Check if the Python object is a wrapped Epetra_IntVector
  int res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_EIV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    result =
      reinterpret_cast< Teuchos::RCP< Epetra_IntVector > * >(argp);
    return result;
  }

  //
  // Get the default communicator. N.B.: it is important that this
  // getComm() method be called after the code just above. If the
  // python script is shutting down and destroying objects, this
  // function could get called by the python destructor for an
  // Epetra.IntVector object. In this case, the default Teuchos
  // communicator might already be destroyed. If that is the case,
  // though, then the code just above will do the conversion and
  // return before getting here.
  const Teuchos::RCP< const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

#ifdef HAVE_PYTRILINOS_DOMI
  //
  // Check if the Python object is a wrapped Domi::MDVector<int>
  *newmem = 0;
  res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_DMDV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    dmdv_rcp =
      *reinterpret_cast< Teuchos::RCP< Domi::MDVector<int> > * >(argp);
    try
    {
      smartresult = dmdv_rcp->getEpetraIntVectorView();
      *newmem = *newmem | SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_IntVector >(smartresult);
    return result;
  }
  //
  // Check if the Python object supports the DistArray Protocol
  if (PyObject_HasAttrString(pyobj, "__distarray__"))
  {
    try
    {
      if (!(distarray = PyObject_CallMethod(pyobj, (char*) "__distarray__", (char*) "")))
        return NULL;
      DistArrayProtocol dap(distarray);
      dmdv_rcp = convertToMDVector<int>(comm, dap);
      *newmem = SWIG_CAST_NEW_MEMORY;
      Py_DECREF(distarray);
    }
    catch (PythonException & e)
    {
      e.restore();
      return NULL;
    }
    try
    {
      smartresult = dmdv_rcp->getEpetraIntVectorView();
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_IntVector >(smartresult);
    return result;
  }
#endif

  //
  // Check if the environment is serial, and if so, check if the
  // Python object is a NumPy array
  if (comm->getSize() == 1)
  {
    if (PyArray_Check(pyobj))
    {
      Epetra_SerialComm comm = Epetra_SerialComm();
      PyArrayObject * array =
        (PyArrayObject*) PyArray_ContiguousFromObject(pyobj, NPY_INT, 0, 0);
      if (!array) return NULL;
      const int totalLength = PyArray_Size((PyObject*)array);
      int * data = (int*) PyArray_DATA(array);
      Epetra_Map map(totalLength, 0, comm);
      smartresult = Teuchos::rcp(new Epetra_IntVector(Copy, map, data));
      *newmem = SWIG_CAST_NEW_MEMORY;
      result = new Teuchos::RCP< Epetra_IntVector >(smartresult);
      return result;
    }
  }
  //
  // If we get to this point, then none of our known converters will
  // work, so it is time to throw an exception.
  PyErr_Format(PyExc_TypeError, "Could not convert argument of type '%s'\n"
               "to an Epetra_IntVector",
               convertPyStringToChar(PyObject_Str(PyObject_Type(pyobj))));
  throw PythonException();
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_MultiVector > *
convertPythonToEpetraMultiVector(PyObject * pyobj,
                                 int * newmem)
{
  // SWIG initialization
  static swig_type_info * swig_EMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector >*");
  static swig_type_info * swig_DMDV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Domi::MDVector<double> >*");
  //
  // Result objects
  void *argp = 0;
  PyObject * distarray = 0;
  Teuchos::RCP< Epetra_MultiVector > smartresult;
  Teuchos::RCP< Epetra_MultiVector > * result;
#ifdef HAVE_PYTRILINOS_DOMI
  Teuchos::RCP< Domi::MDVector<double> > dmdv_rcp;
#endif
  *newmem = 0;
  //
  // Check if the Python object is a wrapped Epetra_MultiVector
  int res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_EMV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    result =
      reinterpret_cast< Teuchos::RCP< Epetra_MultiVector > * >(argp);
    return result;
  }

  //
  // Get the default communicator. N.B.: it is important that this
  // getComm() method be called after the code just above. If the
  // python script is shutting down and destroying objects, this
  // function could get called by the python destructor for an
  // Epetra.MultiVector object. In this case, the default Teuchos
  // communicator might already be destroyed. If that is the case,
  // though, then the code just above will do the conversion and
  // return before getting here.
  const Teuchos::RCP< const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

#ifdef HAVE_PYTRILINOS_DOMI
  //
  // Check if the Python object is a wrapped Domi::MDVector<double>
  *newmem = 0;
  res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_DMDV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    dmdv_rcp =
      *reinterpret_cast< Teuchos::RCP< Domi::MDVector<double> > * >(argp);
    try
    {
      smartresult = dmdv_rcp->getEpetraMultiVectorView();
      *newmem = *newmem | SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_MultiVector >(smartresult);
    return result;
  }
  //
  // Check if the Python object supports the DistArray Protocol
  if (PyObject_HasAttrString(pyobj, "__distarray__"))
  {
    try
    {
      if (!(distarray = PyObject_CallMethod(pyobj, (char*) "__distarray__", (char*) "")))
        return NULL;
      DistArrayProtocol dap(distarray);
      dmdv_rcp = convertToMDVector<double>(comm, dap);
      *newmem = SWIG_CAST_NEW_MEMORY;
      Py_DECREF(distarray);
    }
    catch (PythonException & e)
    {
      e.restore();
      return NULL;
    }
    try
    {
      smartresult = dmdv_rcp->getEpetraMultiVectorView();
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_MultiVector >(smartresult);
    return result;
  }
#endif

  //
  // Check if the environment is serial, and if so, check if the
  // Python object is a NumPy array
  if (comm->getSize() == 1)
  {
    if (PyArray_Check(pyobj))
    {
      Epetra_SerialComm comm = Epetra_SerialComm();
      PyArrayObject * array =
        (PyArrayObject*) PyArray_ContiguousFromObject(pyobj, NPY_DOUBLE, 0, 0);
      if (!array) return NULL;
      int numVec, vecLen;
      int ndim = PyArray_NDIM(array);
      if (ndim == 1)
      {
        numVec = 1;
        vecLen = PyArray_DIM(array, 0);
      }
      else
      {
        numVec = PyArray_DIM(array, 0);
        vecLen = 1;
        for (int i=1; i < ndim; ++i) vecLen *= PyArray_DIM(array, i);
      }
      double * data = (double*) PyArray_DATA(array);
      Epetra_Map map(vecLen, 0, comm);
      smartresult =
        Teuchos::rcp(new Epetra_MultiVector(Copy, map, data, vecLen, numVec));
      result = new Teuchos::RCP< Epetra_MultiVector >(smartresult);
      *newmem = SWIG_CAST_NEW_MEMORY;
      return result;
    }
  }
  //
  // If we get to this point, then none of our known converters will
  // work, so it is time to set a Python error
  PyErr_Format(PyExc_TypeError, "Could not convert argument of type '%s'\n"
               "to an Epetra_MultiVector",
               convertPyStringToChar(PyObject_Str(PyObject_Type(pyobj))));
  return NULL;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Vector > *
convertPythonToEpetraVector(PyObject * pyobj,
                            int * newmem)
{
  // SWIG initialization
  static swig_type_info * swig_EV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector >*");
  static swig_type_info * swig_DMDV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Domi::MDVector<double> >*");
  //
  // Result objects
  void *argp = 0;
  PyObject * distarray = 0;
  Teuchos::RCP< Epetra_Vector > smartresult;
  Teuchos::RCP< Epetra_Vector > * result;
#ifdef HAVE_PYTRILINOS_DOMI
  Teuchos::RCP< Domi::MDVector<double> > dmdv_rcp;
#endif
  *newmem = 0;
  //
  // Check if the Python object is a wrapped Epetra_Vector
  int res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_EV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    result =
      reinterpret_cast< Teuchos::RCP< Epetra_Vector > * >(argp);
    return result;
  }

  //
  // Get the default communicator. N.B.: it is important that this
  // getComm() method be called after the code just above. If the
  // python script is shutting down and destroying objects, this
  // function could get called by the python destructor for an
  // Epetra.Vector object. In this case, the default Teuchos
  // communicator might already be destroyed. If that is the case,
  // though, then the code just above will do the conversion and
  // return before getting here.
  const Teuchos::RCP< const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

#ifdef HAVE_PYTRILINOS_DOMI
  //
  // Check if the Python object is a wrapped Domi::MDVector<double>
  *newmem = 0;
  res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_DMDV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    dmdv_rcp =
      *reinterpret_cast< Teuchos::RCP< Domi::MDVector<double> > * >(argp);
    try
    {
      smartresult = dmdv_rcp->getEpetraVectorView();
      *newmem = *newmem | SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_Vector >(smartresult);
    return result;
  }
  //
  // Check if the Python object supports the DistArray Protocol
  if (PyObject_HasAttrString(pyobj, "__distarray__"))
  {
    try
    {
      if (!(distarray = PyObject_CallMethod(pyobj, (char*) "__distarray__", (char*) "")))
        return NULL;
      DistArrayProtocol dap(distarray);
      dmdv_rcp = convertToMDVector<double>(comm, dap);
      *newmem = SWIG_CAST_NEW_MEMORY;
      Py_DECREF(distarray);
    }
    catch (PythonException & e)
    {
      e.restore();
      return NULL;
    }
    try
    {
      smartresult = dmdv_rcp->getEpetraVectorView();
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Epetra_Vector >(smartresult);
    return result;
  }
#endif

  //
  // Check if the environment is serial, and if so, check if the
  // Python object is a NumPy array
  if (comm->getSize() == 1)
  {
    if (PyArray_Check(pyobj))
    {
      Epetra_SerialComm comm = Epetra_SerialComm();
      PyArrayObject * array =
        (PyArrayObject*) PyArray_ContiguousFromObject(pyobj, NPY_DOUBLE, 0, 0);
      if (!array) return NULL;
      const int totalLength = PyArray_Size((PyObject*)array);
      double * data = (double*) PyArray_DATA(array);
      Epetra_Map map(totalLength, 0, comm);
      smartresult = Teuchos::rcp(new Epetra_Vector(Copy, map, data));
      result = new Teuchos::RCP< Epetra_Vector >(smartresult);
      *newmem = SWIG_CAST_NEW_MEMORY;
      return result;
    }
  }
  //
  // If we get to this point, then none of our known converters will
  // work, so it is time to throw an exception.
  PyErr_Format(PyExc_TypeError, "Could not convert argument of type '%s'\n"
               "to an Epetra_Vector",
               convertPyStringToChar(PyObject_Str(PyObject_Type(pyobj))));
  throw PythonException();
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< Epetra_Operator > *eo)
{
  // SWIG initialization
  static swig_type_info *swig_EO_ptr   =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator        >*");
  // static swig_type_info *swig_EFCO_ptr =
  //   SWIG_TypeQuery("Teuchos::RCP< Epetra_FastCrsOperator >*");
  static swig_type_info *swig_EIO_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_InvOperator     >*");
  static swig_type_info *swig_ERM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_RowMatrix       >*");
  static swig_type_info *swig_EBRM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_BasicRowMatrix  >*");
  static swig_type_info *swig_ECM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_CrsMatrix       >*");
  // static swig_type_info *swig_EMM_ptr  =
  //   SWIG_TypeQuery("Teuchos::RCP< Epetra_MsrMatrix       >*");
  static swig_type_info *swig_EVM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrMatrix       >*");
  static swig_type_info *swig_EVRM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrRowMatrix    >*");
  static swig_type_info *swig_EFVM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_FEVbrMatrix     >*");
  static swig_type_info *swig_EFCM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_FECrsMatrix     >*");
  static swig_type_info *swig_EJM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_JadMatrix       >*");
  //
  // Attempt to downcast to Epetra_VbrRowMatrix
  Teuchos::RCP< Epetra_VbrRowMatrix > *evrm = new
    Teuchos::RCP< Epetra_VbrRowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_VbrRowMatrix >(*eo));
  if (evrm->is_null()) delete evrm;
  else return SWIG_NewPointerObj((void*)evrm, swig_EVRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FEVbrMatrix
  Teuchos::RCP< Epetra_FEVbrMatrix > *efvm = new
    Teuchos::RCP< Epetra_FEVbrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_FEVbrMatrix >(*eo));
  if (efvm->is_null()) delete efvm;
  else return SWIG_NewPointerObj((void*)efvm, swig_EFVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FECrsMatrix
  Teuchos::RCP< Epetra_FECrsMatrix > *efcm = new
    Teuchos::RCP< Epetra_FECrsMatrix >(Teuchos::rcp_dynamic_cast< Epetra_FECrsMatrix >(*eo));
  if (efcm->is_null()) delete efcm;
  else return SWIG_NewPointerObj((void*)efcm, swig_EFCM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_JadMatrix
  Teuchos::RCP< Epetra_JadMatrix > *ejm = new
    Teuchos::RCP< Epetra_JadMatrix >(Teuchos::rcp_dynamic_cast< Epetra_JadMatrix >(*eo));
  if (ejm->is_null()) delete ejm;
  else return SWIG_NewPointerObj((void*)ejm, swig_EJM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_BasicRowMatrix
  Teuchos::RCP< Epetra_BasicRowMatrix > *ebrm = new
    Teuchos::RCP< Epetra_BasicRowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_BasicRowMatrix >(*eo));
  if (ebrm->is_null()) delete ebrm;
  else return SWIG_NewPointerObj((void*)ebrm, swig_EBRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_CrsMatrix
  Teuchos::RCP< Epetra_CrsMatrix > *ecm = new
    Teuchos::RCP< Epetra_CrsMatrix >(Teuchos::rcp_dynamic_cast< Epetra_CrsMatrix >(*eo));
  if (ecm->is_null()) delete ecm;
  else return SWIG_NewPointerObj((void*)ecm, swig_ECM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_MsrMatrix
  // Teuchos::RCP< Epetra_MsrMatrix > *emm = new
  //   Teuchos::RCP< Epetra_MsrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_MsrMatrix >(*eo));
  // if (emm->is_null()) delete emm;
  // else return SWIG_NewPointerObj((void*)emm, swig_EMM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< Epetra_VbrMatrix > *evm = new
    Teuchos::RCP< Epetra_VbrMatrix >(Teuchos::rcp_dynamic_cast< Epetra_VbrMatrix >(*eo));
  if (evm->is_null()) delete evm;
  else return SWIG_NewPointerObj((void*)evm, swig_EVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_RowMatrix
  Teuchos::RCP< Epetra_RowMatrix > *erm = new
    Teuchos::RCP< Epetra_RowMatrix >(Teuchos::rcp_dynamic_cast< Epetra_RowMatrix >(*eo));
  if (erm->is_null()) delete erm;
  else return SWIG_NewPointerObj((void*)erm, swig_ERM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_InvOperator
  Teuchos::RCP< Epetra_InvOperator > *eio = new
    Teuchos::RCP< Epetra_InvOperator >(Teuchos::rcp_dynamic_cast< Epetra_InvOperator >(*eo));
  if (eio->is_null()) delete eio;
  else return SWIG_NewPointerObj((void*)eio, swig_EIO_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FastCrsOperator
  // Teuchos::RCP< Epetra_FastCrsOperator > *efco = new
  //   Teuchos::RCP< Epetra_FastCrsOperator >(Teuchos::rcp_dynamic_cast< Epetra_FastCrsOperator >(*eo));
  // if (efco->is_null()) delete efco;
  // else return SWIG_NewPointerObj((void*)efco, swig_EFCO_ptr, SWIG_POINTER_OWN);
  //
  // If downcast is not possible, return python Epetra_Operator
  // **NOTE** Always make a copy!
  Teuchos::RCP< Epetra_Operator > *eocopy = new Teuchos::RCP< Epetra_Operator >(*eo);
  return SWIG_NewPointerObj((void*)eocopy, swig_EO_ptr, SWIG_POINTER_OWN);
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertEpetraOperatorToPython(const Teuchos::RCP< const Epetra_Operator > *ceo)
{
  // SWIG initialization
  static swig_type_info *swig_EO_ptr   =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator        >*");
  // static swig_type_info *swig_EFCO_ptr =
  //   SWIG_TypeQuery("Teuchos::RCP< Epetra_FastCrsOperator >*");
  static swig_type_info *swig_EIO_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_InvOperator     >*");
  static swig_type_info *swig_ERM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_RowMatrix       >*");
  static swig_type_info *swig_EBRM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_BasicRowMatrix  >*");
  static swig_type_info *swig_ECM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_CrsMatrix       >*");
  // static swig_type_info *swig_EMM_ptr  =
  //   SWIG_TypeQuery("Teuchos::RCP< Epetra_MsrMatrix       >*");
  static swig_type_info *swig_EVM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrMatrix       >*");
  static swig_type_info *swig_EVRM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_VbrRowMatrix    >*");
  static swig_type_info *swig_EFVM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_FEVbrMatrix     >*");
  static swig_type_info *swig_EFCM_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_FECrsMatrix     >*");
  static swig_type_info *swig_EJM_ptr  =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_JadMatrix       >*");
  //
  // Cast const-ness away
  Teuchos::RCP< Epetra_Operator > *eo = new
    Teuchos::RCP< Epetra_Operator >(Teuchos::rcp_const_cast< Epetra_Operator >(*ceo));
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< const Epetra_VbrRowMatrix > *evrm = new
    Teuchos::RCP< const Epetra_VbrRowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_VbrRowMatrix >(*eo));
  if (evrm->is_null()) delete evrm;
  else return SWIG_NewPointerObj((void*)evrm, swig_EVRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FEVbrMatrix
  Teuchos::RCP< const Epetra_FEVbrMatrix > *efvm = new
    Teuchos::RCP< const Epetra_FEVbrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_FEVbrMatrix >(*eo));
  if (efvm->is_null()) delete efvm;
  else return SWIG_NewPointerObj((void*)efvm, swig_EFVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FECrsMatrix
  Teuchos::RCP< const Epetra_FECrsMatrix > *efcm = new
    Teuchos::RCP< const Epetra_FECrsMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_FECrsMatrix >(*eo));
  if (efcm->is_null()) delete efcm;
  else return SWIG_NewPointerObj((void*)efcm, swig_EFCM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_JadMatrix
  Teuchos::RCP< const Epetra_JadMatrix > *ejm = new
    Teuchos::RCP< const Epetra_JadMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_JadMatrix >(*eo));
  if (ejm->is_null()) delete ejm;
  else return SWIG_NewPointerObj((void*)ejm, swig_EJM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_BasicRowMatrix
  Teuchos::RCP< const Epetra_BasicRowMatrix > *ebrm = new
    Teuchos::RCP< const Epetra_BasicRowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_BasicRowMatrix >(*eo));
  if (ebrm->is_null()) delete ebrm;
  else return SWIG_NewPointerObj((void*)ebrm, swig_EBRM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_CrsMatrix
  Teuchos::RCP< const Epetra_CrsMatrix > *ecm = new
    Teuchos::RCP< const Epetra_CrsMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_CrsMatrix >(*eo));
  if (ecm->is_null()) delete ecm;
  else return SWIG_NewPointerObj((void*)ecm, swig_ECM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_MsrMatrix
  // Teuchos::RCP< const Epetra_MsrMatrix > *emm = new
  //   Teuchos::RCP< const Epetra_MsrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_MsrMatrix >(*eo));
  // if (emm->is_null()) delete emm;
  // else return SWIG_NewPointerObj((void*)emm, swig_EMM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_VbrMatrix
  Teuchos::RCP< const Epetra_VbrMatrix > *evm = new
    Teuchos::RCP< const Epetra_VbrMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_VbrMatrix >(*eo));
  if (evm->is_null()) delete evm;
  else return SWIG_NewPointerObj((void*)evm, swig_EVM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_RowMatrix
  Teuchos::RCP< const Epetra_RowMatrix > *erm = new
    Teuchos::RCP< const Epetra_RowMatrix >(Teuchos::rcp_dynamic_cast< const Epetra_RowMatrix >(*eo));
  if (erm->is_null()) delete erm;
  else return SWIG_NewPointerObj((void*)erm, swig_ERM_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_InvOperator
  Teuchos::RCP< const Epetra_InvOperator > *eio = new
    Teuchos::RCP< const Epetra_InvOperator >(Teuchos::rcp_dynamic_cast< const Epetra_InvOperator >(*eo));
  if (eio->is_null()) delete eio;
  else return SWIG_NewPointerObj((void*)eio, swig_EIO_ptr, SWIG_POINTER_OWN);
  //
  // Attempt to downcast to Epetra_FastCrsOperator
  // Teuchos::RCP< const Epetra_FastCrsOperator > *efco = new
  //   Teuchos::RCP< const Epetra_FastCrsOperator >(Teuchos::rcp_dynamic_cast< const Epetra_FastCrsOperator >(*eo));
  // if (efco->is_null()) delete efco;
  // else return SWIG_NewPointerObj((void*)efco, swig_EFCO_ptr, SWIG_POINTER_OWN);
  //
  // If downcast is not possible, return python Epetra_Operator
  return SWIG_NewPointerObj((void*)eo, swig_EO_ptr, SWIG_POINTER_OWN);
}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Vector >
getEpetraVectorObjectAttr(PyObject   * object,
                          CONST char * name)
{
  static swig_type_info * swig_EV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value,
                                                    &argp,
                                                    swig_EV_ptr,
                                                    0,
                                                    &newmem)))
  {
    PyErr_Format(PyExc_TypeError,
                 "Attribute '%s' is not of type Epetra.Vector",
                 name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP< Epetra_Vector > result =
    *reinterpret_cast< Teuchos::RCP< Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_Vector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Epetra_Vector >
getConstEpetraVectorObjectAttr(PyObject   * object,
                               CONST char * name)
{
  static swig_type_info * swig_EV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError,
                 "Attribute '%s' is not of type Epetra.Vector",
                 name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP< const Epetra_Vector > result =
    *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Epetra_Vector >
getConstEpetraVectorItemObjectAttr(PyObject   * object,
                                   CONST char * name,
                                   int          i)
{
  static swig_type_info * swig_EV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Vector > *");
  void * argp;
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(item, &argp, swig_EV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError,
                 "Attribute '%s' is not tuple of type Epetra.Vector",
                 name);
    Py_DECREF(item);
    throw PythonException();
  }
  Teuchos::RCP< const Epetra_Vector > result =
    *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(argp);
  Py_DECREF(item);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_MultiVector >
getEpetraMultiVectorObjectAttr(PyObject   * object,
                               CONST char * name)
{
  static swig_type_info * swig_EMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_MultiVector > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value, &argp, swig_EMV_ptr, 0, &newmem)))
  {
    PyErr_Format(PyExc_TypeError,
                 "Attribute '%s' is not of type Epetra.MultiVector",
                 name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_MultiVector > result =
    *reinterpret_cast< Teuchos::RCP< Epetra_MultiVector > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_MultiVector > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Operator >
getEpetraOperatorObjectAttr(PyObject   * object,
                            CONST char * name)
{
  static swig_type_info * swig_EO_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Epetra_Operator > *");
  void * argp;
  PyObject * value = PyObject_GetAttrString(object, name);
  int newmem = 0;
  if (!SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value,
                                                    &argp,
                                                    swig_EO_ptr,
                                                    0,
                                                    &newmem)))
  {
    PyErr_Format(PyExc_TypeError,
                 "Attribute '%s' is not of type Epetra.Operator",
                 name);
    Py_DECREF(value);
    throw PythonException();
  }
  Teuchos::RCP<Epetra_Operator > result =
    *reinterpret_cast< Teuchos::RCP< Epetra_Operator > * >(argp);
  if (newmem)
    delete reinterpret_cast< Teuchos::RCP< Epetra_Operator > * >(argp);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////

PyObject *
convertToDimData(const Epetra_BlockMap & ebm,
                 int   extraDim)
{
  // Initialization
  PyObject * dim_data    = NULL;
  PyObject * dim_dict    = NULL;
  PyObject * indices     = NULL;
  npy_intp   dims        = 1;
  int        currentDim  = 0;

  // Get the Epetra_Comm
  const Epetra_Comm & comm = ebm.Comm();

  // Get the number of dimensions.  The vector data constitutes one
  // dimension.  If the extraDim is greater than one, then that
  // constitutes a second possible dimension.  If there is a constant
  // element size and that size is greater than one, then that
  // constitutes a third possible dimension.  If the element size is
  // variable, throw an exception, as PyTrilinos and the DistArray
  // Protocol do not handle that case.
  Py_ssize_t ndim     = 1;
  long       elemSize = 1;
  if (extraDim > 1) ndim += 1;
  if (ebm.ConstantElementSize())
  {
    if (ebm.ElementSize() > 1)
    {
      ndim += 1;
      elemSize = ebm.ElementSize();
    }
  }
  else
  {
    PyErr_SetString(PyExc_ValueError,
                    "The given Epetra_BLockMap has variable element sizes, "
                    "which is not currently supported by PyTrilinos for the "
                    "DistArray Protocol.");
    goto fail;
  }

  // Allocate the dim_data return object, which is a tuple of length
  // ndim
  dim_data = PyTuple_New(ndim);
  if (!dim_data) goto fail;

  // If we have an extra dimension argument greater than one, then
  // define a dimension associated with the multiple vectors
  if (extraDim > 1)
  {
    dim_dict = PyDict_New();
    if (!dim_dict) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("n")) == -1) goto fail;

    if (PyDict_SetItemString(dim_dict,
                             "size",
                             PyInt_FromLong(extraDim)) == -1) goto fail;
    // This function steals a reference from dim_dict
    PyTuple_SET_ITEM(dim_data, currentDim, dim_dict);
    currentDim += 1;
  }

  // If we have a constant element size greater than one, then define
  // a dimension associated with the element size
  if (elemSize > 1)
  {
    dim_dict = PyDict_New();
    if (!dim_dict) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("n")) == -1) goto fail;

    if (PyDict_SetItemString(dim_dict,
                             "size",
                             PyInt_FromLong(elemSize)) == -1) goto fail;
    // This function steals a reference from dim_dict
    PyTuple_SET_ITEM(dim_data, currentDim, dim_dict);
    currentDim += 1;
  }

  // Allocate the dimension data dictionary for the distributed points
  // and assign the common key values
  dim_dict = PyDict_New();
  if (!dim_dict) goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "size",
                           PyInt_FromLong(ebm.NumMyElements())) == -1) goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "proc_grid_size",
                           PyInt_FromLong(comm.NumProc())) == -1) goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "proc_grid_rank",
                           PyInt_FromLong(comm.MyPID())) == -1) goto fail;

  // Determine the type of the dimension data, either block "b" or
  // unstructured "u", set the "dist_type" key and other keys required
  // according to dist_type.
  if (ebm.LinearMap())
  {
    // LinearMap == true corresponds to DistArray Protocol dist_type
    // == "b" (block)
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("b")) == -1) goto fail;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
    long long minMyGID = (long long) ebm.MinMyGID();
    long long maxMyGID = (long long) ebm.MaxMyGID();
#else
    long long minMyGID = ebm.MinMyGID64();
    long long maxMyGID = ebm.MaxMyGID64();
#endif
    if (PyDict_SetItemString(dim_dict,
                             "start",
                             PyInt_FromLong(minMyGID)) == -1) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "stop",
                             PyInt_FromLong(maxMyGID+1)) == -1) goto fail;
  }
  else
  {
    // LinearMap == false corresponds to DistArray Protocol dist_type
    // == "u" (unstructured)
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("u")) == -1) goto fail;
    dims    = ebm.NumMyElements();
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
    indices = PyArray_SimpleNewFromData(1,
                                        &dims,
                                        NPY_LONG,
                                        (void*)ebm.MyGlobalElements64());
#else
    indices = PyArray_SimpleNewFromData(1,
                                        &dims,
                                        NPY_INT,
                                        (void*)ebm.MyGlobalElements());
#endif
    if (!indices) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "indices",
                             indices) == -1) goto fail;
    Py_DECREF(indices);
    indices = NULL;
  }

  // Put the dimension dictionary into the dim_data tuple
  PyTuple_SET_ITEM(dim_data, currentDim, dim_dict);

  // Return the dim_data tuple
  return dim_data;

  // Error handling
  fail:
  Py_XDECREF(dim_data);
  Py_XDECREF(dim_dict);
  Py_XDECREF(indices);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertToDistArray(const Epetra_IntVector & eiv)
{
  // Initialization
  PyObject   * dap      = NULL;
  PyObject   * dim_data = NULL;
  PyObject   * dim_dict = NULL;
  PyObject   * size     = NULL;
  PyObject   * buffer   = NULL;
  Py_ssize_t   ndim     = 1;
  npy_intp     dims[3];

  // Get the underlying Epetra_BlockMap
  const Epetra_BlockMap & ebm = eiv.Map();

  // Allocate the DistArray object and set the version key value
  dap = PyDict_New();
  if (!dap) goto fail;
  if (PyDict_SetItemString(dap,
                           "__version__",
                           PyString_FromString("0.9.0")) == -1) goto fail;

  // Get the Dimension Data and the number of dimensions.  If the
  // underlying Epetra_BlockMap has variable element sizes, an error
  // will be detected here.
  dim_data = convertToDimData(ebm);
  if (!dim_data) goto fail;
  ndim = PyTuple_Size(dim_data);

  // Assign the Dimension Data key value.
  if (PyDict_SetItemString(dap,
                           "dim_data",
                           dim_data) == -1) goto fail;

  // Extract the buffer dimensions from the Dimension Data, construct
  // the buffer and assign the buffer key value
  for (Py_ssize_t i = 0; i < ndim; ++i)
  {
    dim_dict = PyTuple_GetItem(dim_data, i);
    if (!dim_dict) goto fail;
    size = PyDict_GetItemString(dim_dict, "size");
    if (!size) goto fail;
    dims[i] = PyInt_AsLong(size);
    if (PyErr_Occurred()) goto fail;
  }
  buffer = PyArray_SimpleNewFromData(ndim,
                                     dims,
                                     NPY_INT,
                                     (void*)eiv.Values());
  if (!buffer) goto fail;
  if (PyDict_SetItemString(dap,
                           "buffer",
                           buffer) == -1) goto fail;

  // Return the DistArray Protocol object
  return dap;

  // Error handling
  fail:
  Py_XDECREF(dap);
  Py_XDECREF(dim_data);
  Py_XDECREF(dim_dict);
  Py_XDECREF(buffer);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject *
convertToDistArray(const Epetra_MultiVector & emv)
{
  // Initialization
  PyObject   * dap      = NULL;
  PyObject   * dim_data = NULL;
  PyObject   * dim_dict = NULL;
  PyObject   * size     = NULL;
  PyObject   * buffer   = NULL;
  Py_ssize_t   ndim     = 1;
  npy_intp     dims[3];

  // Get the underlying Epetra_BlockMap
  const Epetra_BlockMap & ebm = emv.Map();

  // Allocate the DistArray object and set the version key value
  dap = PyDict_New();
  if (!dap) goto fail;
  if (PyDict_SetItemString(dap,
                           "__version__",
                           PyString_FromString("0.9.0")) == -1) goto fail;

  // Get the Dimension Data and the number of dimensions.  If the
  // underlying Epetra_BlockMap has variable element sizes, an error
  // will be detected here.
  dim_data = convertToDimData(ebm,
                              emv.NumVectors());
  if (!dim_data) goto fail;
  ndim = PyTuple_Size(dim_data);

  // Assign the Dimension Data key value.
  if (PyDict_SetItemString(dap,
                           "dim_data",
                           dim_data) == -1) goto fail;

  // Extract the buffer dimensions from the Dimension Data, construct
  // the buffer and assign the buffer key value
  for (Py_ssize_t i = 0; i < ndim; ++i)
  {
    dim_dict = PyTuple_GetItem(dim_data, i);
    if (!dim_dict) goto fail;
    size = PyDict_GetItemString(dim_dict, "size");
    if (!size) goto fail;
    dims[i] = PyInt_AsLong(size);
    if (PyErr_Occurred()) goto fail;
  }
  buffer = PyArray_SimpleNewFromData(ndim,
                                     dims,
                                     NPY_DOUBLE,
                                     (void*)emv[0]);
  if (!buffer) goto fail;
  if (PyDict_SetItemString(dap,
                           "buffer",
                           buffer) == -1) goto fail;

  // Return the DistArray Protocol object
  return dap;

  // Error handling
  fail:
  Py_XDECREF(dap);
  Py_XDECREF(dim_data);
  Py_XDECREF(dim_dict);
  Py_XDECREF(buffer);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

}  // Namespace PyTrilinos
