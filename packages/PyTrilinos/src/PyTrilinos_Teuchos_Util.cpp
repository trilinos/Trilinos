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

// Include files
#include "Python3Compat.hpp"
#include "PyTrilinos_Teuchos_Util.hpp"
#include "swigpyrun.h"

// Backward compatibility for python < 2.5
#if (PY_VERSION_HEX < 0x02050000)
typedef int Py_ssize_t;
#endif

// ****************************************************************** //

namespace PyTrilinos
{

// ****************************************************************** //

bool setPythonParameter(Teuchos::ParameterList & plist,
			const std::string      & name,
			PyObject               * value)
{
  static swig_type_info * swig_TPL_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Teuchos::ParameterList >*");
  void * argp;
  int newmem = 0;

  // Boolean values
  if (PyBool_Check(value))
  {
    if (value == Py_True) plist.set(name,true );
    else                  plist.set(name,false);
  }

  // Integer values
  else if (PyInt_Check(value))
  {
    plist.set(name, (int)PyInt_AsLong(value));
  }

  // Floating point values
  else if (PyFloat_Check(value))
  {
    plist.set(name, PyFloat_AsDouble(value));
  }

  // Unicode values
  else if (PyUnicode_Check(value))
  {
    PyObject * pyBytes = PyUnicode_AsASCIIString(value);
    if (!pyBytes) return false;
    plist.set(name, std::string(PyBytes_AsString(pyBytes)));
    Py_DECREF(pyBytes);
  }

  // String values
  else if (PyString_Check(value))
  {
    plist.set(name, std::string(convertPyStringToChar(value)));
  }

  // None object not allowed: this is a python type not usable by
  // Trilinos solver packages, so we reserve it for the
  // getPythonParameter() function to indicate that the requested
  // parameter does not exist in the given Teuchos::ParameterList.
  // For logic reasons, this check must come before the check for
  // Teuchos::ParameterList
  else if (value == Py_None)
  {
    return false;
  }

  // Dictionary values.  This must come before the check for Python
  // sequences, because the Python ditionary is a sequence.
  else if (PyDict_Check(value))
  {  
    // Convert the python dictionary to a Teuchos::ParameterList
    Teuchos::ParameterList * sublist = pyDictToNewParameterList(value);

    // Store the Teuchos::ParameterList
    plist.set(name,*sublist);
    delete sublist;
  }

  // Teuchos::ParameterList values
  else if (SWIG_CheckState(SWIG_Python_ConvertPtrAndOwn(value,
                                                        &argp,
                                                        swig_TPL_ptr,
                                                        0,
                                                        &newmem)))
  {
    if (newmem & SWIG_CAST_NEW_MEMORY)
    { 
      Teuchos::RCP< Teuchos::ParameterList > tempshared =
	*reinterpret_cast< Teuchos::RCP< Teuchos::ParameterList > * >(argp);
      delete reinterpret_cast< Teuchos::RCP< Teuchos::ParameterList > * >(argp);
      plist.set(name, *(tempshared.get()));
    }
    else
    {
      Teuchos::RCP< Teuchos::ParameterList > * smartarg =
	reinterpret_cast< Teuchos::RCP< Teuchos::ParameterList > * >(argp);
      if (smartarg) plist.set(name, *(smartarg->get()));
    }
  }

  // NumPy arrays and non-dictionary Python sequences
  else if (PyArray_Check(value) || PySequence_Check(value))
  {
    PyObject * pyArray =
      PyArray_CheckFromAny(value,
                           NULL,
                           1,
                           1,
                           NPY_ARRAY_DEFAULT | NPY_ARRAY_NOTSWAPPED,
                           NULL);
    if (!pyArray) return false;
    // if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_BOOL)
    // {
    //   Teuchos::Array< bool > tArray;
    //   copyNumPyToTeuchosArray(pyArray, tArray);
    //   plist.set(name, tArray);
    // }
    // else if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_INT)
    if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_INT)
    {
      Teuchos::Array< int > tArray;
      copyNumPyToTeuchosArray(pyArray, tArray);
      plist.set(name, tArray);
    }
    else if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_LONG)
    {
      Teuchos::Array< long > tArray;
      copyNumPyToTeuchosArray(pyArray, tArray);
      plist.set(name, tArray);
    }
    else if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_FLOAT)
    {
      Teuchos::Array< float > tArray;
      copyNumPyToTeuchosArray(pyArray, tArray);
      plist.set(name, tArray);
    }
    else if (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_DOUBLE)
    {
      Teuchos::Array< double > tArray;
      copyNumPyToTeuchosArray(pyArray, tArray);
      plist.set(name, tArray);
    }
    else if ((PyArray_TYPE((PyArrayObject*) pyArray) == NPY_STRING) ||
             (PyArray_TYPE((PyArrayObject*) pyArray) == NPY_UNICODE))
    {
      Teuchos::Array< std::string > tArray;
      copyNumPyToTeuchosArray(pyArray, tArray);
      plist.set(name, tArray);
    }
    else
    {
      // Unsupported data type
      if (pyArray != value) Py_DECREF(pyArray);
      return false;
    }
  }

  // All other value types are unsupported
  else
  {
    return false;
  }

  // Successful type conversion
  return true;
}    // setPythonParameter

// **************************************************************** //

PyObject * getPythonParameter(const Teuchos::ParameterList & plist,
			      const std::string            & name)
{
  static swig_type_info * swig_TPL_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Teuchos::ParameterList >*");

  // If parameter does not exist, return None
  if (!plist.isParameter(name)) return Py_BuildValue("");

  // Get the parameter entry.  I now deal with the Teuchos::ParameterEntry
  // objects so that I can query the Teuchos::ParameterList without setting
  // the "used" flag to true.
  const Teuchos::ParameterEntry * entry = plist.getEntryPtr(name);

  // Boolean parameter values
  if (entry->isType< bool >())
  {
    bool value = Teuchos::any_cast< bool >(entry->getAny(false));
    return PyBool_FromLong((long)value);
  }
  // Integer parameter values
  else if (entry->isType< int >())
  {
    int value = Teuchos::any_cast< int >(entry->getAny(false));
    return PyInt_FromLong((long)value);
  }
  // Double parameter values
  else if (entry->isType< double >())
  {
    double value = Teuchos::any_cast< double >(entry->getAny(false));
    return PyFloat_FromDouble(value);
  }
  // String parameter values
  else if (entry->isType< std::string >())
  {
    std::string value = Teuchos::any_cast< std::string >(entry->getAny(false));
    return PyString_FromString(value.c_str());
  }
  // Char * parameter values
  else if (entry->isType< char * >())
  {
    char * value = Teuchos::any_cast< char * >(entry->getAny(false));
    return PyString_FromString(value);
  }
  // Teuchos::ParameterList values
  else if (entry->isList())
  {
    Teuchos::ParameterList value;
    Teuchos::RCP< Teuchos::ParameterList > * valuercp =
      new Teuchos::RCP< Teuchos::ParameterList >(
        new Teuchos::ParameterList(Teuchos::any_cast< Teuchos::ParameterList >(entry->getAny(false))));
    return SWIG_NewPointerObj((void*)valuercp, swig_TPL_ptr, SWIG_POINTER_OWN);
  }
  // Teuchos::Array values
  else if (entry->isArray())
  {
    // try
    // {
    //   Teuchos::Array< bool > tArray =
    //     Teuchos::any_cast< Teuchos::Array< bool > >(entry->getAny(false));
    //   return copyTeuchosArrayToNumPy(tArray);
    // }
    // catch(Teuchos::bad_any_cast &e)
    // {
      try
      {
        Teuchos::Array< int > tArray =
          Teuchos::any_cast< Teuchos::Array< int > >(entry->getAny(false));
        return copyTeuchosArrayToNumPy(tArray);
      }
      catch(Teuchos::bad_any_cast &e)
      {
        try
        {
          Teuchos::Array< long > tArray =
            Teuchos::any_cast< Teuchos::Array< long > >(entry->getAny(false));
          return copyTeuchosArrayToNumPy(tArray);
        }
        catch(Teuchos::bad_any_cast &e)
        {
          try
          {
            Teuchos::Array< float > tArray =
              Teuchos::any_cast< Teuchos::Array< float > >(entry->getAny(false));
            return copyTeuchosArrayToNumPy(tArray);
          }
          catch(Teuchos::bad_any_cast &e)
          {
            try
            {
              Teuchos::Array< double > tArray =
                Teuchos::any_cast< Teuchos::Array< double > >(entry->getAny(false));
              return copyTeuchosArrayToNumPy(tArray);
            }
            catch(Teuchos::bad_any_cast &e)
            {
              try
              {
                Teuchos::Array< std::string > tArray =
                  Teuchos::any_cast< Teuchos::Array< std::string > >(entry->getAny(false));
                return copyTeuchosArrayToNumPy(tArray);
              }
              catch(Teuchos::bad_any_cast &e)
              {
                // Teuchos::Arrays of type other than int or double are
                // currently unsupported
                return NULL;
              }
            }
          }
        }
      }
    // }
  }

  // All  other types are unsupported
  return NULL;
}    // getPythonParameter

// **************************************************************** //

bool isEquivalent(PyObject                     * dict,
                  const Teuchos::ParameterList & plist)
{
  PyObject  * key   = NULL;
  PyObject  * value = NULL;
  PyObject  * param = NULL;
  Py_ssize_t  pos   = 0;
  std::string name;

  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict)) goto fail;

  // Check that all entries in Teuchos::ParameterList are also in the
  // python dictionary
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++)
  {
    name  = plist.name(i);
    value = PyDict_GetItemString(dict,name.c_str());
    if (value == NULL) goto fail;
    if (plist.isSublist(name))
    {
      if (!isEquivalent(value, plist.sublist(name)))
	goto fail;
    }
    else
    {
      param = getPythonParameter(plist,name.c_str());
      if (param == NULL) goto fail;
      if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
      Py_DECREF(param);
    }
  }
  // Check that all entries in the python dictionary are also in
  // the Teuchos::ParameterList
  while (PyDict_Next(dict, &pos, &key, &value))
  {
    if (!PyString_Check(key)) goto fail;
    name = std::string(convertPyStringToChar(key));
    if (!plist.isParameter(name)) goto fail;
    if (plist.isSublist(name))
    {
      if (!isEquivalent(value, plist.sublist(name))) goto fail;
    }
    else
    {
      param = getPythonParameter(plist, name);
      if (param == NULL) goto fail;
      if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
      Py_DECREF(param);
    }
  }
  // All checks passed
  return true;

  fail:
  Py_XDECREF(param);
  return false;
}    // isEquivalent

// **************************************************************** //

bool updatePyDictWithParameterList(PyObject                     * dict,
				   const Teuchos::ParameterList & plist,
				   ResponseToIllegalParameters    flag)
{
  static char  illegalParam[ ] = "Illegal Parameters";
  PyObject   * value   = NULL;
  PyObject   * param   = NULL;
  bool         result  = true;
  const char * nameStr = NULL;
  std::string  name;

  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
    goto fail;
  }

  // Iterate over all entries in Teuchos::ParameterList and ensure they are
  // mirrored in the python dictionary
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++)
  {
    name    = plist.name(i);
    nameStr = name.c_str();
    param   = getPythonParameter(plist,nameStr);
    value   = PyDict_GetItemString(dict,nameStr);

    // If param is NULL, then behavior is determined by flag
    if (param == NULL)
    {
      switch (flag)
      {
      case ignore:
	break;
      case storeNames:
	result = false;
	value = PyDict_GetItemString(dict, illegalParam);
	if (value == NULL)
	{
	  PyDict_SetItemString(dict, illegalParam, Py_BuildValue("(s)", nameStr));
	}
	else
	{
	  if (!PyTuple_Check(value))
	  {
	    PyErr_Format(PyExc_ValueError, "Parameter '%s' has unexpected type",
			 illegalParam);
	    goto fail;
	  }
	  PyDict_SetItemString(dict, illegalParam,
			       PySequence_Concat(value, Py_BuildValue("(s)", nameStr)));
	}
	break;
      case raiseError:
	PyErr_Format(PyExc_ValueError, "Parameter '%s' is of unsupported type", nameStr);
	goto fail;
      default:
	PyErr_Format(PyExc_RuntimeError, "Unexpected enumeration encountered");
	goto fail;
      }
    }
    else
    {
      // If param is a sublist, mirror with a dictionary by calling
      // this routine recursively
      if (plist.isSublist(name))
      {
	if (value == NULL) value = PyDict_New();
	else if (!PyDict_Check(value))
	{
	  Py_DECREF(value);
	  value = PyDict_New();
	}
	result = result and 
	  updatePyDictWithParameterList(value, plist.sublist(name));
	PyDict_SetItemString(dict,nameStr,value);
      }

      // Else synchronize the dictionary value to the parameter
      else PyDict_SetItemString(dict,nameStr,param);
    }
    Py_XDECREF(param);
  }
  return result;
  fail:
  return false;
}    // updatePyDictWithParameterList

// **************************************************************** //

bool updateParameterListWithPyDict(PyObject                  * dict,
				   Teuchos::ParameterList    & plist,
				   ResponseToIllegalParameters flag)
{
  static char illegalKey[ ] = "Illegal Keys";
  PyObject  * key    = NULL;
  PyObject  * value  = NULL;
  Py_ssize_t  pos    = 0;
  bool        result = true;
  std::string name;

  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
    goto fail;
  }

  // Iterate over all items in the python dictionary and ensure they
  // are synchronized with the Teuchos::ParameterList
  while (PyDict_Next(dict, &pos, &key, &value))
  {

    // If the key is not a string, we can't synchronize
    if (!PyUnicode_Check(key) && !PyString_Check(key))
    {
      PyErr_SetString(PyExc_TypeError, "Encountered non-string key in dictionary");
      goto fail;
    }

    if (PyUnicode_Check(key))
    {
      PyObject * pyBytes = PyUnicode_AsASCIIString(key);
      if (!pyBytes) goto fail;
      name = std::string(PyBytes_AsString(pyBytes));
      Py_DECREF(pyBytes);
    }
    else
    {
      name = std::string(convertPyStringToChar(key));
    }
    if (!setPythonParameter(plist, name, value))
    {
      // If value is not settable, behavior is determined by flag
      switch (flag)
      {
      case ignore:
	break;
      case storeNames:
	result = false;
	if (plist.isParameter(illegalKey))
	{
	  const Teuchos::ParameterEntry * entry = plist.getEntryPtr(name);
	  if (entry->isType< std::string >())
	  {
	    std::string names = Teuchos::any_cast< std::string >(entry->getAny(false));
	    plist.set(illegalKey, names + std::string(", ") + name);
	  }
	  else
	  {
	    PyErr_Format(PyExc_TypeError, "Parameter '%s' has unexpected type",
			 illegalKey);
	    goto fail;
	  }
	}
	else
	{
	  plist.set(illegalKey, name);
	}
	break;
      case raiseError:
	PyErr_Format(PyExc_ValueError, "Parameter '%s' has unsupported value",
		     illegalKey);
	goto fail;
      default:
	PyErr_Format(PyExc_RuntimeError, "Unexpected enumeration encountered");
	goto fail;
      }
    }
  }
  return result;
  fail:
  return false;
}    // updateParameterListWithPyDict

// **************************************************************** //

bool synchronizeParameters(PyObject                  * dict,
			   Teuchos::ParameterList    & plist,
			   ResponseToIllegalParameters flag)
{
  bool result = updatePyDictWithParameterList(dict,plist,flag);
  result      = result && updateParameterListWithPyDict(dict,plist,flag);
  return result;
}    // synchronizeParameters

// **************************************************************** //

Teuchos::ParameterList *
pyDictToNewParameterList(PyObject                  * dict,
			 ResponseToIllegalParameters flag)
{
  Teuchos::ParameterList * plist = 0;
  // The dict pointer must point to a dictionary
  if (!PyDict_Check(dict))
  {
    PyErr_SetString(PyExc_ValueError, "Expected a dictionary");
    goto fail;
  }

  // Create a new Teuchos::ParameterList and synchronize it with the python
  // dictionary
  plist = new Teuchos::ParameterList();
  if (!updateParameterListWithPyDict(dict,*plist,flag))
  {
    // If update failed, behavior is determined by flag
    switch (flag)
    {
    case ignore:
    case storeNames:
      break;
    case raiseError:
    default:
      delete plist;
      goto fail;
    }
  }
  return plist;
  fail:
  return NULL;
}    // pyDictToNewParameterList

// **************************************************************** //

PyObject * parameterListToNewPyDict(const Teuchos::ParameterList & plist,
				    ResponseToIllegalParameters    flag)
{

  // Create a new dictionary and synchronize it with the Teuchos::ParameterList
  PyObject * dict = PyDict_New();
  if (!updatePyDictWithParameterList(dict,plist))
  {

    // If update failed, behavior is determined by flag
    switch (flag)
    {
    case ignore:
    case storeNames:
      break;
    case raiseError:
    default:
      Py_XDECREF(dict);
      goto fail;
    }
  }
  return dict;
  fail:
  return NULL;
}    // parameterListToNewPyDict

}    // Namespace PyTrilinos
