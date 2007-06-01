
#include "python_ParameterEntry.hpp"
#include "Teuchos_FILEstream.hpp"
#include <boost/python.hpp>

using namespace boost::python;
using namespace Teuchos;

/******************************************************************/
// Type method: return python type of requested parameter.  Replaces
// templated C++ isType() methods.
PyObject * type( ParameterList& self,const std::string & name) {
  PyObject * value = python_plist_tools::getPythonParameter(self,name);
  // Type not supported
  if (value == NULL) {
PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
goto fail;
  }
  // Name not found
  else if (value == Py_None) {
PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
goto fail;
  }
  // Name found and type supported
  return PyObject_Type(value);
fail:
  return NULL;
}

//////////////////////////////////////////////
// ** The following methods are added to ** //
// ** give ParameterList a PyDict "feel" ** //
//////////////////////////////////////////////

/******************************************************************/
// Comparison operators
int __cmp__(ParameterList& self, PyObject * obj) {
  PyObject * dict = python_plist_tools::parameterListToNewPyDict(self,ignore);
  int result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_Compare(dict,obj);
  Py_DECREF(dict);
  return result;
fail:
  Py_XDECREF(dict);
  return -2;
}

int __cmp__(ParameterList& self,const ParameterList & ParameterList) {
  PyObject * dict1 = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * dict2 = python_plist_tools::parameterListToNewPyDict(ParameterList,ignore);
  int result = 0;
  if (dict1 == NULL) goto fail;
  if (dict2 == NULL) goto fail;
  result = PyObject_Compare(dict1,dict2);
  Py_DECREF(dict1);
  Py_DECREF(dict2);
  return result;
fail:
  Py_XDECREF(dict1);
  Py_XDECREF(dict2);
  return -2;
}

/******************************************************************/
// Contains operator
int __contains__(ParameterList& self,const std::string & name) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * keys   = 0;
  PyObject * keyStr = 0;
  int result        = 0;
  if (dict == NULL) goto fail;
  // There is a simpler form in python 2.4, but this should be
  // backward-compatible
  keys   = PyDict_Keys(dict);
  keyStr = PyString_FromString(name.c_str());
  result = PySequence_Contains(keys,keyStr);
  Py_DECREF(dict  );
  Py_DECREF(keys  );
  Py_DECREF(keyStr);
  return result;
fail:
  Py_XDECREF(dict);
  return result;
}

/******************************************************************/
// Equals operators
PyObject * __eq__(ParameterList& self, PyObject * obj) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_RichCompare(dict,obj,Py_EQ);
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

PyObject * __eq__(ParameterList& self, const ParameterList & ParameterList) {
  PyObject * dict1  = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * dict2  = python_plist_tools::parameterListToNewPyDict(ParameterList,ignore);
  PyObject * result = 0;
  if (dict1 == NULL) goto fail;
  if (dict2 == NULL) goto fail;
  result = PyObject_RichCompare(dict1,dict2,Py_EQ);
  Py_DECREF(dict1);
  Py_DECREF(dict2);
  return result;
fail:
  Py_XDECREF(dict1);
  Py_XDECREF(dict2);
  return NULL;
}

/******************************************************************/
/******************************************************************/
// __iter__ method
PyObject * __iter__(ParameterList& self) {
  PyObject * dict = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * iter = 0;
  if (dict == NULL) goto fail;
  iter = PyObject_GetIter(PyDict_Keys(dict));
  Py_DECREF(dict);
  return iter;
fail:
  return NULL;
}

/******************************************************************/
// Length operator
int __len__(ParameterList& self) {
  PyObject * dict = python_plist_tools::parameterListToNewPyDict(self,ignore);
  int len = 0;
  if (dict == NULL) goto fail;
  len = PyDict_Size(dict);
  Py_DECREF(dict);
  return len;
fail:
  return -1;
}

/******************************************************************/
// Not equals operators
PyObject * __ne__(ParameterList& self, PyObject * obj) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_RichCompare(dict,obj,Py_NE);
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

PyObject * __ne__(ParameterList& self, const ParameterList & ParameterList) {
  PyObject * dict1  = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * dict2  = python_plist_tools::parameterListToNewPyDict(ParameterList,ignore);
  PyObject * result = 0;
  if (dict1 == NULL) goto fail;
  if (dict2 == NULL) goto fail;
  result = PyObject_RichCompare(dict1,dict2,Py_NE);
  Py_DECREF(dict1);
  Py_DECREF(dict2);
  return result;
fail:
  Py_XDECREF(dict1);
  Py_XDECREF(dict2);
  return NULL;
}

/******************************************************************/
/******************************************************************/
// String representation method
PyObject * __repr__(ParameterList& self) {
  std::string reprStr;
  PyObject * dict    = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * dictStr = 0;
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  dictStr = PyObject_Str(dict);
  reprStr = std::string("ParameterList(") + 
        std::string(PyString_AsString(dictStr)) +
        std::string(")");
  result = PyString_FromString(reprStr.c_str());
  Py_DECREF(dict   );
  Py_DECREF(dictStr);
  return result;
fail:
  Py_XDECREF(dict);
  return NULL;
}

/******************************************************************/
// String conversion method
PyObject * __str__(ParameterList& self) {
  PyObject * dict = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * str  = 0;
  if (dict == NULL) goto fail;
  str = PyObject_Str(dict);
  Py_DECREF(dict);
  return str;
fail:
  return NULL;
}

/******************************************************************/
// Has_key method
int has_key(ParameterList& self, const std::string & name) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * keys   = 0;
  PyObject * keyStr = 0;
  int result        = 0;
  if (dict == NULL) goto fail;
  // There is a simpler form in python 2.4, but this should be
  // backward-compatible
  keys   = PyDict_Keys(dict);
  keyStr = PyString_FromString(name.c_str());
  result = PySequence_Contains(keys,keyStr);
  Py_DECREF(dict  );
  Py_DECREF(keys  );
  Py_DECREF(keyStr);
  return result;
fail:
  return -2;
}

/******************************************************************/
// Items method
PyObject * items(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyDict_Items(dict);
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// Iteritems method
PyObject * iteritems(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_GetIter(PyDict_Items(dict));
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// Iterkeys method
PyObject * iterkeys(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_GetIter(PyDict_Keys(dict));
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// Itervalues method
PyObject * itervalues(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyObject_GetIter(PyDict_Values(dict));
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// Keys method
PyObject * keys(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyDict_Keys(dict);
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// Update methods
void update(ParameterList& self, PyObject * dict, bool strict=true) {
  ResponseToIllegalParameters flag;
  if (strict) flag = raiseError;
  else        flag = storeNames;
  python_plist_tools::updateParameterListWithPyDict(dict,self,flag);
}

void update(ParameterList& self,const ParameterList & ParameterList) {
  self.setParameters(ParameterList);
}

/******************************************************************/
// Values method
PyObject * values(ParameterList& self) {
  PyObject * dict   = python_plist_tools::parameterListToNewPyDict(self,ignore);
  PyObject * result = 0;
  if (dict == NULL) goto fail;
  result = PyDict_Values(dict);
  Py_DECREF(dict);
  return result;
fail:
  return NULL;
}

/******************************************************************/
// AsDict method: return a dictionary equivalent to the
// ParameterList
PyObject * asDict(ParameterList& self) {
  return python_plist_tools::parameterListToNewPyDict(self,storeNames);
}

/******************************************************************/

PyObject * _myprint(ParameterList& self,PyObject * pf, int indent, bool showTypes,
	      bool showFlags) {
  PyObject * returnObject = pf;
  // No arguments
  if (pf==NULL) {
self.print(std::cout,indent,showTypes,showFlags);
returnObject = Py_None;
  }
  // Given non-file pf argument
  else {
if (!PyFile_Check(pf)) {
  PyErr_SetString(PyExc_IOError, "_print() method expects a file object");
  goto fail;
}
// Given file pf argument
else {
  std::FILE *f = PyFile_AsFile(pf);
  Teuchos::FILEstream buffer(f);
  std::ostream os(&buffer);
  self.print(os,indent,showTypes,showFlags);
}
  }
  Py_INCREF(returnObject);
  return returnObject;
fail:
  return NULL;
}
