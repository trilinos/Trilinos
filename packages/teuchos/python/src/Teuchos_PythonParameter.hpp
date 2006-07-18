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
      PyDictParameterList * pdplist = NULL;
      // Convert the python dictionary to a PyDictParameterList
      PyDictParameterList sublist = PyDictParameterList(value,name);
      if (PyErr_Occurred()) return false;
      // Try to cast plist to PyDictParameterList; if you can, store the sublist
      pdplist = dynamic_cast<PyDictParameterList *>(&plist);
      if (pdplist) pdplist->set(name, sublist);
      // Cast failed -> plist is a ParameterList, so store an upcast sublist
      else plist.sublist(name) = dynamic_cast<ParameterList &>(sublist);
    }

    // None object
    else if (value == Py_None) {
      return false;
    }

    // PyDictParameterList values
    else if (SWIG_CheckState(SWIG_Python_ConvertPtr(value, &argp, swig_TPDPL_ptr, 0))) {
      PyDictParameterList *arg = reinterpret_cast<PyDictParameterList *>(argp);
      plist.set(name, *arg);
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

    static swig_type_info * swig_TPL_ptr   = SWIG_TypeQuery("Teuchos::ParameterList *"      );
    static swig_type_info * swig_TPDPL_ptr = SWIG_TypeQuery("Teuchos::PyDictParameterList *");

    // Check for parameter existence
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

    // Char* parameter values
    else if (isParameterType<char*>(plist,name)) {
      char* value = getParameter<char*>(plist,name);
      return PyString_FromString(value);
    }

    // PyDictParameterList values
    else if (isParameterType<PyDictParameterList>(plist,name)) {
      //PyDictParameterList * value = getParameter<PyDictParameterList>(plist, name);
      //const PyDictParameterList * value = plist.getPtr<const PyDictParameterList>(name);
      const PyDictParameterList & value = plist.get<const PyDictParameterList>(name);
      return SWIG_NewPointerObj((void*) &value, swig_TPDPL_ptr, 0);
    }

    // ParameterList values
    else if (isParameterType<ParameterList>(plist,name)) {
      //ParameterList * value = getParameter<ParameterList>(plist,name);
      //const ParameterList * value = plist.getPtr<const ParameterList>(name);
      const ParameterList & value = plist.get<const ParameterList>(name);
      return SWIG_NewPointerObj((void*) &value, swig_TPL_ptr, 0);
    }

    // Unsupported type
    return NULL;

  }    // getPythonParameter

}    // namespace Teuchos

#endif TEUCHOS_PYTHONPARAMETER
