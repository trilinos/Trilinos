#ifndef TEUCHOS_PYTHONPARAMETER
#define TEUCHOS_PYTHONPARAMETER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_PyDictParameterList.hpp"

namespace Teuchos {

  // Teuchos::ParameterList can support parameters of any type, but the
  // python wrappers need a subset of types supported a-priori.  Since
  // this subset shows up in many places, the setPythonParameter()
  // function is provided here that handles conversion from PyObjects to
  // the supported types, assigning it to the specified ParameterList.
  // It returns true if the conversion was made, false if the requested
  // conversion is not supported.
  bool setPythonParameter(ParameterList     & plist,
			  const std::string & name,
			  PyObject          * value) {

    static swig_type_info * swig_TPL_ptr   = SWIG_TypeQuery("Teuchos::ParameterList *"      );
    static swig_type_info * swig_TPDPL_ptr = SWIG_TypeQuery("Teuchos::PyDictParameterList *");
    void * argp;

    // Boolean values
    if (PyBool_Check(value)) {
      if (value == Py_True) plist.set(name,true );
      else                  plist.set(name,false);
    }

    // Integer values
    else if (PyInt_Check(value)) {
      plist.set(name, (int)PyInt_AsLong(value));
    }

    // Floating point values
    else if (PyFloat_Check(value)) {
      plist.set(name, PyFloat_AsDouble(value));
    }

    // String values
    else if (PyString_Check(value)) {
      plist.set(name, PyString_AsString(value));
    }

    // Dictionary values
    else if (PyDict_Check(value)) {

      // Convert the python dictionary to a PyDictParameterList
      const PyDictParameterList pdplist(value,name);
      if (PyErr_Occurred()) return false;

      // Convert the PyDictParameterList to a regular ParameterList
      const ParameterList sublist(pdplist);

      // Store the ParameterList
      plist.set(name,sublist);
    }

    // None object
    else if (value == Py_None) {
      return false;
    }

    // PyDictParameterList values get stored as ParameterLists
    else if (SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_TPDPL_ptr, 0))) {
      PyDictParameterList *arg = reinterpret_cast<PyDictParameterList *>(argp);
      ParameterList sublist(*arg);
      plist.set(name, sublist);
    }

    // ParameterList values
    else if (SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_TPL_ptr, 0))) {
      ParameterList *arg = reinterpret_cast<ParameterList *>(argp);
      plist.set(name, *arg);
    }

    // Unsupported value types
    else {
      return false;
    }

    // Successful type conversion
    return true;
  }

  // The getPythonParameter() function is the get counterpart to
  // setPythonParameter().  If the requested parameter name does not
  // exist, None is returned (a type that is guaranteed not to be
  // supported).  If the name exists and its type is supported, it is
  // returned as a python object.  If the name exists, but the type is
  // not supported, NULL is returned, to indicate an error.
  PyObject * getPythonParameter(const ParameterList & plist,
				const std::string   & name) {

    static swig_type_info * swig_TPL_ptr = SWIG_TypeQuery("Teuchos::ParameterList *");

    // If parameter does not exist, return None
    if (!plist.isParameter(name)) return Py_BuildValue("");

    // Boolean parameter values
    if (isParameterType<bool>(plist,name)) {
      bool value = getParameter<bool>(plist,name);
      return PyBool_FromLong((long)value);
    }

    // Integer parameter values
    else if (isParameterType<int>(plist,name)) {
      int value = getParameter<int>(plist,name);
      return PyInt_FromLong((long)value);
    }

    // Double parameter values
    else if (isParameterType<double>(plist,name)) {
      double value = getParameter<double>(plist,name);
      return PyFloat_FromDouble(value);
    }

    // String parameter values
    else if (isParameterType<std::string>(plist,name)) {
      std::string value = getParameter<std::string>(plist,name);
      return PyString_FromString(value.c_str());
    }

    // Char * parameter values
    else if (isParameterType<char *>(plist,name)) {
      char * value = getParameter<char *>(plist,name);
      return PyString_FromString(value);
    }

    // ParameterList values
    else if (isParameterType<ParameterList>(plist,name)) {
      const ParameterList & value = plist.sublist(name);
      return SWIG_NewPointerObj((void*) &value, swig_TPL_ptr, 0);
    }

    // Unsupported type
    return NULL;

  }    // getPythonParameter

  // Function isEquivalent() is a utility for determining whether a
  // python dictionary is functionally equivalent to a ParameterList.
  // It supports interpreting ParameterList sublists as nested python
  // dictionaries, so it calls itself recursively.
  bool isEquivalent(PyObject * dict, const ParameterList & plist) {
    PyObject * key   = NULL;
    PyObject * value = NULL;
    PyObject * param = NULL;
    int        pos   = 0;
    string     name;

    // The dict pointer must point to a dictionary
    if (!PyDict_Check(dict)) return false;

    // Check that all entries in ParameterList are also in the
    // python dictionary
    for (ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
      name = plist.name(i);
      value = PyDict_GetItemString(dict ,name.c_str());
      if (value == NULL) goto fail;
      if (plist.isSublist(name)) {
	if (!isEquivalent(value, plist.sublist(name)))
	  goto fail;
      }
      else {
	param = getPythonParameter(  plist,name.c_str());
	if (param == NULL) goto fail;
	if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
	Py_DECREF(param);
      }
    }
    // Check that all entries in the python dictionary are also in
    // the ParameterList
    while (PyDict_Next(dict, &pos, &key, &value)) {
      if (!PyString_Check(key)) goto fail;
      name  = string(PyString_AsString(key));
      if (!plist.isParameter(name)) goto fail;
      if (plist.isSublist(name)) {
	if (!isEquivalent(value, plist.sublist(name)))
	  goto fail;
      }
      else {
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

}    // namespace Teuchos

#endif TEUCHOS_PYTHONPARAMETER
