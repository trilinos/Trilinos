#ifndef PYTHON_PARAMETER_ENTRY
#define PYTHON_PARAMETER_ENTRY

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_RefCountPtr.hpp"
// #include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_ParameterList.hpp"


  // The getPythonParameter() function is the get counterpart to
  // setPythonParameter().  It takes a ParameterList and string name
  // as input.  If the requested parameter name does not exist, a new
  // reference to None is returned (a type that is guaranteed not to
  // be supported by setPythonParameter()).  If the name exists and
  // its type is supported, it is returned as a new reference to a
  // python object.  If the name exists, but the type is not
  // supported, NULL is returned, to indicate an error.  All returned
  // python object pointers are new references.  This function is
  // coded in such a way that the ParameterList "used" flags are not
  // altered.


// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;

// ****************************************************************** //

// The following enumeration is used as an optional flag in certain
// routines to indicate how the routine is supposed to react to
// illegal parameters.

enum ResponseToIllegalParameters {
                raiseError,
			    ignore,
			    storeNames };

// ****************************************************************** //

#if (PY_VERSION_HEX < 0x02050000)
typedef int Py_ssize_t;
#endif


class python_plist_tools{
public:
    static PyObject * getPythonParameter(const Teuchos::ParameterList & plist,
  				const std::string   & name) {
        
      // If parameter does not exist, return None
      if (!plist.isParameter(name)) return Py_BuildValue("");

      // Get the parameter entry.  I now deal with the ParameterEntry
      // objects so that I can query the ParameterList without setting
      // the used flag to true.
      const Teuchos::ParameterEntry * entry = plist.getEntryPtr(name);

      // Boolean parameter values
      if (entry->isType<bool>()) {
        bool value = Teuchos::any_cast<bool>(entry->getAny(false));
        return PyBool_FromLong((long)value);
      }

      // Integer parameter values
      else if (entry->isType<int>()) {
        int value = Teuchos::any_cast<int>(entry->getAny(false));
        return PyInt_FromLong((long)value);
      }

      // Double parameter values
      else if (entry->isType<double>()) {
        double value = Teuchos::any_cast<double>(entry->getAny(false));
        return PyFloat_FromDouble(value);
      }

      // String parameter values
      else if (entry->isType<std::string>()) {
        std::string value = Teuchos::any_cast<std::string>(entry->getAny(false));
        return PyString_FromString(value.c_str());
      }

      // Char * parameter values
      else if (entry->isType<char *>()) {
        char * value = Teuchos::any_cast<char *>(entry->getAny(false));
        return PyString_FromString(value);
      }

      // ParameterList values
      else if (entry->isList()) {
          
        const Teuchos::ParameterList & value = Teuchos::getValue<Teuchos::ParameterList>(*entry);
        object o( &value );
        return o.ptr();
        //return SWIG_NewPointerObj((void*) &value, swig_TPL_ptr, 0);
      }

      // All  other types are unsupported
      PyErr_SetString( PyExc_TypeError, "unsupported type for python plist" );
      throw_error_already_set();
      return NULL;

    }    // getPythonParameter

      
      // ########################################################################################
      // ########################################################################################
      
      static bool setPythonParameter(Teuchos::ParameterList     & plist,
    			  const std::string & name,
    			  PyObject*           value) {

        // handle<> o(borrowed(value) );
        // object ob( o );
        // Teuchos::ParameterList *pp = dynamic_cast< Teuchos::ParameterList >( ob.ptr() );
        
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
          plist.set(name, std::string(PyString_AsString(value)));
        }

        // Dictionary values
        else if (PyDict_Check(value)) {
          // Convert the python dictionary to a ParameterList
          Teuchos::ParameterList * sublist = pyDictToNewParameterList(value);

          // Store the ParameterList
          plist.set(name,*sublist);
          delete sublist;
        }

        // None object not allowed: this is a python type not usable by
        // Trilinos solver packages, so we reserve it for the
        // getPythonParameter() function to indicate that the requested
        // parameter does not exist in the given ParameterList
        else if (value == Py_None) {
            PyErr_SetString( PyExc_TypeError, "unsupported type for python plist" );
            throw_error_already_set();
            return false;
        }

        // ParameterList values
        // else if ( pp != NULL ) {
        //   plist.set(name, *pp);
        // }

        // All other value types are unsupported
        else {
            PyErr_SetString( PyExc_TypeError, "unsupported type for python plist" );
            throw_error_already_set();
            return false;
        }

        // Successful type conversion
        return true;
      }    // setPythonParameter
      
      // ########################################################################################
      // ########################################################################################
      
        static bool isEquivalent(PyObject * dict, const Teuchos::ParameterList & plist) {
          PyObject * key   = NULL;
          PyObject * value = NULL;
          PyObject * param = NULL;
          Py_ssize_t pos   = 0;
          string     name;

          // The dict pointer must point to a dictionary
          if (!PyDict_Check(dict)) goto fail;

          // Check that all entries in ParameterList are also in the
          // python dictionary
          for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
            name  = plist.name(i);
            value = PyDict_GetItemString(dict,name.c_str());
            if (value == NULL) goto fail;
            if (plist.isSublist(name)) {
      	if (!isEquivalent(value, plist.sublist(name)))
      	  goto fail;
            }
            else {
      	param = getPythonParameter( plist,name.c_str() ) ;
      	if (param == NULL) goto fail;
      	if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
      	Py_DECREF(param);
            }
          }
          // Check that all entries in the python dictionary are also in
          // the ParameterList
          while (PyDict_Next(dict, &pos, &key, &value)) {
            if (!PyString_Check(key)) goto fail;
            name = string(PyString_AsString(key));
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

      // ########################################################################################
      // ########################################################################################

        // **************************************************************** //

        // Function updatePyDictWithParameterList() takes all of the entries
        // in a ParameterList and updates the given python dictionary to
        // reflect the same values.  If the given python object is not a
        // dictionary, or any of the ParameterList entries are of
        // unsupported type, the function returns false.

        static bool updatePyDictWithParameterList(PyObject * dict, const Teuchos::ParameterList & plist,
      				     ResponseToIllegalParameters flag=raiseError) {
          static char  illegalParam[ ] = "Illegal Parameters";
          PyObject   * value   = NULL;
          PyObject   * param   = NULL;
          bool         result  = true;
          const char * nameStr = NULL;
          string       name;

          // The dict pointer must point to a dictionary
          if (!PyDict_Check(dict)) {
            PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
            throw_error_already_set();
            goto fail;
          }

          // Iterate over all entries in ParameterList and ensure they are
          // mirrored in the python dictionary
          for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
            name    = plist.name(i);
            nameStr = name.c_str();
            param   = getPythonParameter(plist,nameStr);
            value   = PyDict_GetItemString(dict,nameStr);

            // If param is NULL, then behavior is determined by flag
            if (param == NULL) {
      	switch (flag) {
      	case ignore:
      	  break;
      	case storeNames:
      	  result = false;
      	  value = PyDict_GetItemString(dict, illegalParam);
      	  if (value == NULL) {
      	    PyDict_SetItemString(dict, illegalParam, Py_BuildValue("(s)", nameStr));
      	  }
      	  else {
      	    if (!PyTuple_Check(value)) {
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
            else {

      	// If param is a sublist, mirror with a dictionary by calling
      	// this routine recursively
      	if (plist.isSublist(name)) {
      	  if (value == NULL) value = PyDict_New();
      	  else if (!PyDict_Check(value)) {
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

        // Function updateParameterListWithPyDict() takes all of the entries
        // in a python dictionary and updates the given ParameterList to
        // reflect the same values.  If the given python object is not a
        // dictionary, or if any of the dictionary keys are not strings, or
        // if any of the dictionary values are of unsupported type, then the
        // function returns false.

        static bool updateParameterListWithPyDict(PyObject * dict, Teuchos::ParameterList & plist,
      				     ResponseToIllegalParameters flag=raiseError) {
          static char illegalKey[ ] = "Illegal Keys";
          PyObject *  key    = NULL;
          PyObject *  value  = NULL;
          Py_ssize_t  pos    = 0;
          bool        result = true;
          string      name;

          // The dict pointer must point to a dictionary
          if (!PyDict_Check(dict)) {
            PyErr_SetString(PyExc_TypeError, "Expected a dictionary");
            throw_error_already_set();
            goto fail;
          }

          // Iterate over all items in the python dictionary and ensure they
          // are synchronized with the ParameterList
          while (PyDict_Next(dict, &pos, &key, &value)) {

            // If the key is not a string, we can't synchronize
            if (!PyString_Check(key)) {
      	PyErr_SetString(PyExc_TypeError, "Encountered non-string key in dictionary");
      	goto fail;
            }

            name = string(PyString_AsString(key));
            if (!setPythonParameter(plist, name, value)) {

      	// If value is not settable, behavior is determined by flag
      	switch (flag) {
      	case ignore:
      	  break;
      	case storeNames:
      	  result = false;
      	  if (plist.isParameter(illegalKey)) {
      	    const Teuchos::ParameterEntry * entry = plist.getEntryPtr(name);
      	    if (entry->isType<string>()) {
      	      string names = Teuchos::any_cast<string>(entry->getAny(false));
      	      plist.set(illegalKey, names + string(", ") + name);
      	    }
      	    else {
      	      PyErr_Format(PyExc_TypeError, "Parameter '%s' has unexpected type",
      			   illegalKey);
      	      goto fail;
      	    }
      	  }
      	  else {
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

        // Function synchronizeParameters() is a function for bringing the
        // given python dictionary and ParameterList into synch with each
        // other.  If a parameter exists for both the ParameterList and the
        // python dictionary, the ParameterList takes precedence.  If the
        // function returns false, it means the given PyObject was not a
        // dictionary or the ParameterList or python dictionary had at least
        // one value of an unsupported type.

        static bool synchronizeParameters(PyObject * dict, Teuchos::ParameterList & plist,
      			     ResponseToIllegalParameters flag=raiseError) {
          bool result = true;
          result = result && updatePyDictWithParameterList(dict,plist,flag);
          result = result && updateParameterListWithPyDict(dict,plist,flag);
          return result;
        }    // synchronizeParameters

        // **************************************************************** //

        // Function pyDictToNewParameterList is a helper function that takes a
        // python dictionary and returns a pointer to an equivalent, new
        // ParameterList.  If dict is not a python dictionary, or dict is
        // not a valid dictionary (non-string keys or unsupported value
        // types) then the function returns NULL.

        static Teuchos::ParameterList * pyDictToNewParameterList(PyObject * dict,
      					   ResponseToIllegalParameters flag=raiseError) {

            Teuchos::ParameterList * plist = 0;
            // The dict pointer must point to a dictionary
            if (!PyDict_Check(dict)) {
                PyErr_SetString(PyExc_ValueError, "Expected a dictionary");
                goto fail;
            }

            // Create a new ParameterList and synchronize it with the python
            // dictionary
            plist = new Teuchos::ParameterList();
            if ( !updateParameterListWithPyDict(dict,*plist,flag) ) {
                // If update failed, behavior is determined by flag
                switch (flag) 
                {
                    case ignore:
                    case storeNames:
                        break;
                    case raiseError:
                        delete plist;
                        goto fail;
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

        // Function parameterListToNewPyDict is a helper function that takes
        // a ParameterList and returns a pointer to an equivalent, new
        // python dictionary.  If the ParameterList contains entries of
        // invalid type, then a python error is raised and NULL is returned.

        static PyObject * parameterListToNewPyDict(const Teuchos::ParameterList & plist,
      				      ResponseToIllegalParameters flag=raiseError) {

            // Create a new dictionary and synchronize it with the ParameterList
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

    

}; // python_plist_tools

#endif


