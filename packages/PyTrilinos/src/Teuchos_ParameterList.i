// -*- c++ -*-

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

%{
// Teuchos includes
#include "Teuchos_any.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_FILEstream.hpp"

// Teuchos python interface includes
#include "Teuchos_PythonParameter.h"
%}

// The python implementation of Teuchos::ParameterList is augmented to
// "play nice" with python dictionaries, and to some extent, behave
// like a dictionary.  A constructor is added that takes a dictionary.
// In fact, any method that takes a ParameterList is extended to also
// take a dictionary.  Most dictionary methods and operators are added
// to the python implementation of ParameterList, so that it can be
// treated as a dictionary.  The exception to this is deletion
// methods, such as pop(), popitem() or __delitem__().  Since
// ParameterList does not support deletion, these methods are not
// implemented.

// Handle the Teuchos::ParameterList:PrintOptions nested class by
// defining it exclusively for SWIG as though it were not nested.
namespace Teuchos
{
class PrintOptions
{
public:
  PrintOptions();
  PrintOptions& indent(int _indent);
  PrintOptions& showTypes(bool _showTypes);
  PrintOptions& showFlags(bool _showFlags);
  PrintOptions& showDoc(bool _showDoc);
  PrintOptions& incrIndent(int indents);
  int indent() const;
  bool showTypes() const;
  bool showFlags() const;
  bool showDoc() const;
  PrintOptions copy() const;
};
%nestedworkaround ParameterList::PrintOptions;
}  // Namespace Teuchos

//Teuchos imports
namespace Teuchos
{
class any;
}
%import "Teuchos_ParameterEntry.hpp"
%import "Teuchos_PythonParameter.h"

/////////////////////////////////////////////////////////////////////
// Override typemaps for ParameterLists to allow PyDicts as input //
/////////////////////////////////////////////////////////////////////
%define %teuchos_rcp_pydict_overrides(CONST, CLASS...)
// Input a plain reference
%typemap(in) CONST CLASS &
(void *argp=0, int res=0, bool cleanup=false, Teuchos::RCP< CONST CLASS > tempshared)
{
  if (PyDict_Check($input))
  {
    $1 = PyTrilinos::pyDictToNewParameterList($input);
    if (!$1) SWIG_fail;
    cleanup = true;
  }
  else
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(Teuchos::RCP< CLASS > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum); 
    }
    if (!argp)
    {
      %argument_nullref("$type", $symname, $argnum);
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      tempshared = *%reinterpret_cast(argp, Teuchos::RCP< CONST CLASS > *);
      delete %reinterpret_cast(argp, Teuchos::RCP< CONST CLASS > *);
      $1 = %const_cast(tempshared.get(), $1_ltype);
    }
    else
    {
      $1 = %const_cast(%reinterpret_cast(argp, Teuchos::RCP< CONST CLASS > *)->get(), $1_ltype);
    }
  }
}
// Perform type checking
%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER,noblock=1)
  CONST CLASS,
  CONST CLASS &,
  CONST CLASS *,
  CONST CLASS *&,
  Teuchos::RCP< CONST CLASS >,
  Teuchos::RCP< CONST CLASS > &,
  Teuchos::RCP< CONST CLASS > *,
  Teuchos::RCP< CONST CLASS > *&
{
  // Accept PyDicts or CLASS instances
  $1 = PyDict_Check($input);
  if (!$1)
  {
    int res = SWIG_ConvertPtr($input, 0, $descriptor(Teuchos::RCP< CLASS > *), 0);
    $1 = SWIG_CheckState(res);
  }
}
// Cleanup
%typemap(freearg) CONST CLASS &
{
  if (cleanup$argnum && $1) delete $1;
}
// Input a Teuchos::RCP by reference
%typemap(in) Teuchos::RCP< CONST CLASS > &
(void *argp, int res = 0, $*1_ltype tempshared)
{
  if (PyDict_Check($input))
  {
    tempshared = Teuchos::rcp(PyTrilinos::pyDictToNewParameterList($input));
    if (tempshared.is_null()) SWIG_fail;
    $1 = &tempshared;
  }
  else
  {
    int newmem = 0;
    res = SWIG_ConvertPtrAndOwn($input, &argp, $descriptor(Teuchos::RCP< CLASS > *),
				%convertptr_flags, &newmem);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum); 
    }
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      if (argp) tempshared = *%reinterpret_cast(argp, $ltype);
      delete %reinterpret_cast(argp, $ltype);
      $1 = &tempshared;
    }
    else
    {
      $1 = (argp) ? %reinterpret_cast(argp, $ltype) : &tempshared;
    }
  }
}
%enddef

////////////////////////////////////
// Teuchos::ParameterList support //
////////////////////////////////////

// Implement internal storage of Teuchos::ParameterList objects via
// Teuchos::RCP<>
%teuchos_rcp(Teuchos::ParameterList)
#define EMPTYHACK
%teuchos_rcp_pydict_overrides(EMPTYHACK, Teuchos::ParameterList)
%teuchos_rcp_pydict_overrides(const,     Teuchos::ParameterList)

%feature("docstring") Teuchos::ParameterList
"The ``ParameterList`` class is an important utility class that is used
by several Trilinos packages for communicating arbitrary-type
parameters between users and packages.

Often, the ``ParameterList`` class is invisible to the user.  It is
analagous to the python dictionary (with the restriction that the
dictionary keys must be strings), and python programmers can provide a
python dictionary wherever a ``ParameterList`` is expected.
``Teuchos`` is imported by other packages that use the
``ParameterList`` class and converts between dictionaries and
``ParameterList`` objects automatically.

The user can create a ``Teuchos.ParameterList`` directly, using the
constructor, ``set`` and ``sublist`` methods, if he so chooses, and
methods that accept ``ParameterList`` objects will work as expected.
It is really just a question of verbosity and elegance that argues in
favor of using a python dictionary.

The python implementation of the ``ParameterList`` class has been
expanded extensively.  Its constructor can accept a python dictionary,
and several methods and operators have been added to the class so that
it behaves somewhat like a dictionary.

C++ ``ParameterList`` objects are designed to support parameters of
arbitrary type.  The python implementation supports a subset of types
*a priori* :

  +-------------------------+-----+-------------------+
  |       Python type       | Dir |     C/C++ type    |
  +-------------------------+-----+-------------------+
  | ``bool``                | <-> | ``bool``          |
  +-------------------------+-----+-------------------+
  | ``int``                 | <-> | ``int``           |
  +-------------------------+-----+-------------------+
  | ``float``               | <-> | ``double``        |
  +-------------------------+-----+-------------------+
  | ``str``                 | <-- | ``char *``        |
  +-------------------------+-----+-------------------+
  | ``str``                 | <-> | ``std::string``   |
  +-------------------------+-----+-------------------+
  | ``dict``                | --> | ``ParameterList`` |
  +-------------------------+-----+-------------------+
  | ``ParameterList``       | <-> | ``ParameterList`` |
  +-------------------------+-----+-------------------+

The C++ ``ParameterList`` class supports ``begin()`` and ``end()``
methods for iterating over the parameters.  These methods are disabled
in the python implementation, in favor of the dictionary iterator
methods: ``__iter__()``, ``iteritems()``, ``iterkeys()`` and
``itervalues()``.

Note that the C++ implementation of the ``ParameterList`` class does
not support parameter deletion.  Therefore, python dictionary methods
that delete items, such as ``pop()`` or ``__delitem__()``, have not
been added to the ``ParameterList`` class."
%feature("docstring")
Teuchos::ParameterList::ParameterList
"If ``dict`` is provided, it must be a dictionary whose keys are all
strings and whose values are all of supported types.  The string name
argument is optional and defaults to ``ANONYMOUS``."
%feature("docstring")
Teuchos::ParameterList::set
"The templated C++ ``set()`` method is replaced in python with a method
that takes a string name and a python object of supported type.  For
example::

    plist = Teuchos.ParameterList()
    plist.set('b',True)
    plist.set('i',10)
    plist.set('f',2.718)
    plist.set('s','Trilinos')
    plist.set('d',{'a':1, 'b':2})"
%feature("docstring")
Teuchos::ParameterList::get
"The templated C++ ``get()`` method is replaced in python with a method
that returns a python object.  For example::

    plist = Teuchos.ParameterList({'f':2.718, 'd':{'a':1. 'b':2}})
    print plist.get('f')
    print plist.get('d')

  will output::

    2.718
    {'a': 1, 'b': 2}"
%feature("docstring")
Teuchos::ParameterList::setParameters
"The ``setParameters()`` method can take either a ``ParameterList`` or
a python dictionary as its argument.  The ``ParameterList`` is updated
to contain all of the entries of the argument."
%feature("docstring")
Teuchos::ParameterList::_print
"Since ``print`` is a python keyword, the ``print()`` C++ method has
been renamed ``_print`` in python.  It takes an optional file argument
that defaults to standard output.  Its output is the same as the C++
implementation."
%feature("autodoc",
"__str__(self) -> string

The ``__str__()`` method returns a string representation of the
``ParameterList`` as though it were a python dictionary.  The python
``eval`` function applied to the output of ``__str__()`` will produce
an equivalent dictionary.")
Teuchos::ParameterList::__str__;
%feature("autodoc",
"__repr__(self) -> string

The ``__repr__()`` method returns the ``__str__()`` output
encapsulated by ``ParameterList(...)``.  The python ``eval`` function
applied to the output of ``__repr__()`` will produce an equivalent
``ParameterList``.")
Teuchos::ParameterList::__repr__;
%feature("autodoc",
"unused(self, pf=None)

The ``unused()`` method in python takes an optional python file object
as its argument, defaulting to standard output.  This specifies where
the output should go.")
Teuchos::ParameterList::unused;
%feature("docstring") Teuchos::ParameterList::type
"Parameter type determination.  With the templated ``isType()`` methods
disabled, type determination of python ``ParameterList`` entries is
accomplished with the ``type()`` method, which returns python type
objects.  For example::

  plist = Teuchos.ParameterList({'b':True, 's':'Trilinos'})
  print plist.type('b')
  print plist.type('s')

results in::

  <type 'bool'>
  <type 'str'>

A non-existent key given as the argument will raise a ``KeyError``
exception."
%feature("docstring")
Teuchos::ParameterList::asDict
"The ``asDict()`` method has been added to ``ParameterList``, which
returns the contents of the ``ParameterList`` converted to a python
dictionary.
"
%feature("docstring")
Teuchos::ParameterList::__eq__
"The ``ParameterList`` equals operator (==)"
%feature("docstring")
Teuchos::ParameterList::__ne__
"The ``ParameterList`` not equals operator (!=)"
%feature("docstring")
Teuchos::ParameterList::__contains__
"The python ``in`` operator works for ``ParameterList`` objects,
searching the parameter names::

  plist = Teuchos.ParameterList({'b':False})
  print 'a' in plist
  print 'b' in plist

produces::

  False
  True"
%feature("docstring")
Teuchos::ParameterList::__len__
"The python ``len()`` function works on ``ParameterList`` objects just
as on python dictionaries::

  plist = Teuchos.ParameterList({'b':True,
                                 'i':10,
                                 'f':2.718,
                                 's':'Trilinos',
                                 'd':{'a':1, 'b':2}})
  print len(plist)

gives::

  5"
%feature("docstring")
Teuchos::ParameterList::__iter__
"To iterate over the parameters in a ``ParameterList``, treat it like a
dictionary::

  plist = Teuchos.ParameterList({'b':True,
                                 'i':10,
                                 'f':2.718,
                                 's':'Trilinos',
                                 'd':{'a':1, 'b':2}})
  for key in plist:
    print key, ':', plist[key]

will result in the output::

  b : True
  d : {'a': 1, 'b': 2}
  f : 2.718
  i : 10
  s : Trilinos

Note that the order of the parameters is somewhat indeterminant, as
with dictionaries, because the iteration object is obtained from an
equivalent dictionary, and dictionaries are ordered by hash
function."
%feature("docstring")
Teuchos::ParameterList::__setitem__
"Like dictionaries, parameters can be set using square brackets::

  plist = Teuchos.ParameterList()
  plist['zero'] = 0"
%feature("docstring")
Teuchos::ParameterList::__getitem__
"Like dictionaries, parameters can be gotten using square brackets::

  plist = Teuchos.ParameterList()
  plist['f'] = 2.718
  e = plist['f']"
%feature("docstring")
Teuchos::ParameterList::update
"An ``update()`` method has been added to the ``ParameterList`` class
that can accept either a ``ParameterList`` or a python dictionary.
Otherwise, it behaves just as the dictionary method, which is
functionally equivalent to the ``setParameters()`` method."
%feature("docstring")
Teuchos::ParameterList::has_key
"Equivalent to the python dictionary has_key() method"
%feature("docstring")
Teuchos::ParameterList::items
"Equivalent to the python dictionary items() method"
%feature("docstring")
Teuchos::ParameterList::iteritems
"Equivalent to the python dictionary iteritems() method"
%feature("docstring")
Teuchos::ParameterList::iterkeys
"Equivalent to the python dictionary iterkeys() method"
%feature("docstring")
Teuchos::ParameterList::itervalues
"Equivalent to the python dictionary itervalues() method"
%feature("docstring")
Teuchos::ParameterList::keys
"Equivalent to the python dictionary keys() method"
%feature("docstring")
Teuchos::ParameterList::values
"Equivalent to the python dictionary values() method"
%extend Teuchos::ParameterList
{
  /******************************************************************/
  // Dictionary constructor
  ParameterList(PyObject * dict, string name = string("ANONYMOUS"))
  {
    Teuchos::ParameterList * plist =
      PyTrilinos::pyDictToNewParameterList(dict, PyTrilinos::raiseError);
    if (plist == NULL) goto fail;

    plist->setName(name);
    return plist;
  fail:
    return NULL;
  }

  /******************************************************************/
  // Set method: accept only python objects as values
  Teuchos::ParameterList & set(const string &name, PyObject *value)
  {
    if (!PyTrilinos::setPythonParameter(*self,name,value))
    {
      PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
      goto fail;
    }
    return *self;
  fail:
    return *self;
  }

  /******************************************************************/
  // Get method: return entries as python objects
  PyObject * get(const string &name, PyObject * default_value=NULL) const
  {
    PyObject * value = PyTrilinos::getPythonParameter(*self, name);
    // Type not supported
    if (value == NULL)
    {
      PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
      goto fail;
    }
    // Name not found
    else if (value == Py_None)
    {
      if (default_value == NULL)
      {
	PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
	goto fail;
      }
      Py_DECREF(value);
      Py_INCREF(default_value);
      return default_value;
    }
    // Type supported and name found
    else return value;
  fail:
    Py_XDECREF(value);
    return NULL;
  }

  /******************************************************************/
  // Print method: change name from "print" to "_print" because
  // "print" is a python keyword
  PyObject * _print(PyObject * pf=NULL, int indent=0, bool showTypes=false,
		    bool showFlags=true)
  {
    PyObject * returnObject = pf;
    // No arguments
    if (pf==NULL)
    {
      self->print(std::cout,indent,showTypes,showFlags);
      returnObject = Py_None;
    }
    // Given non-file pf argument
    else
    {
      if (!PyFile_Check(pf))
      {
	PyErr_SetString(PyExc_IOError, "_print() method expects a file object");
	goto fail;
      }
      // Given file pf argument
      else
      {
	std::FILE *f = PyFile_AsFile(pf);
	Teuchos::FILEstream buffer(f);
	std::ostream os(&buffer);
	self->print(os,indent,showTypes,showFlags);
      }
    }
    Py_INCREF(returnObject);
    return returnObject;
  fail:
    return NULL;
  }

  /******************************************************************/
  // Unused method: take python file objects rather than ostreams
  PyObject * unused(PyObject * pf=NULL)
  {
    // No arguments
    if (pf==NULL)
    {
      self->unused(std::cout);
    }
    // Given non-file pf argument
    else
    {
      if (!PyFile_Check(pf))
      {
	PyErr_SetString(PyExc_IOError, "unused() method expects a file object");
	goto fail;
      }
      // Given file pf argument
      else
      {
	std::FILE *f = PyFile_AsFile(pf);
	Teuchos::FILEstream buffer(f);
	std::ostream os(&buffer);
	self->unused(os);
      }
    }
    return Py_BuildValue("");
  fail:
    return NULL;
  }

  /******************************************************************/
  // Type method: return python type of requested parameter.  Replaces
  // templated C++ isType() methods.
  PyObject * type(const std::string & name)
  {
    PyObject * result = NULL;
    PyObject * value  = PyTrilinos::getPythonParameter(*self,name);
    // Type not supported
    if (value == NULL)
    {
      PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
      goto fail;
    }
    // Name not found
    else if (value == Py_None)
    {
      PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
      goto fail;
    }
    // Name found and type supported
    result = PyObject_Type(value);
    Py_DECREF(value);
    return result;
  fail:
    Py_XDECREF(value);
    return NULL;
  }

  //////////////////////////////////////////////
  // ** The following methods are added to ** //
  // ** give ParameterList a PyDict "feel" ** //
  //////////////////////////////////////////////

  /******************************************************************/
  // Comparison operators
  int __cmp__(PyObject * obj) const
  {
    PyObject * dict = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    int result = 0;
    if (dict == NULL) goto fail;
    result = PyObject_Compare(dict,obj);
    Py_DECREF(dict);
    return result;
  fail:
    Py_XDECREF(dict);
    return -2;
  }

  int __cmp__(const ParameterList & plist) const
  {
    PyObject * dict1 = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    PyObject * dict2 = PyTrilinos::parameterListToNewPyDict(plist,PyTrilinos::ignore);
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
  int __contains__(const std::string & name) const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * __eq__(PyObject * obj) const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    PyObject * result = 0;
    if (dict == NULL) goto fail;
    result = PyObject_RichCompare(dict,obj,Py_EQ);
    Py_DECREF(dict);
    return result;
  fail:
    return NULL;
  }

  PyObject * __eq__(const ParameterList & plist) const
  {
    PyObject * dict1  = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    PyObject * dict2  = PyTrilinos::parameterListToNewPyDict(plist,PyTrilinos::ignore);
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
  // GetItem operator
  PyObject * __getitem__(const std::string & name) const
  {
    // I'm using SWIG's mangling scheme here
    // return Teuchos_ParameterList_get__SWIG_0(self,name);
    return Teuchos_ParameterList_get(self,name);
  }

  /******************************************************************/
  // __iter__ method
  PyObject * __iter__() const
  {
    PyObject * dict = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  int __len__() const
  {
    PyObject * dict = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * __ne__(PyObject * obj) const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    PyObject * result = 0;
    if (dict == NULL) goto fail;
    result = PyObject_RichCompare(dict,obj,Py_NE);
    Py_DECREF(dict);
    return result;
  fail:
    return NULL;
  }

  PyObject * __ne__(const ParameterList & plist) const
  {
    PyObject * dict1  = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
    PyObject * dict2  = PyTrilinos::parameterListToNewPyDict(plist,PyTrilinos::ignore);
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
  // SetItem operator
  void __setitem__(const std::string & name, PyObject * value)
  {
    // I'm using SWIG's mangling scheme here
    Teuchos_ParameterList_set(self,name,value);
  }

  /******************************************************************/
  // String representation method
  PyObject * __repr__() const
  {
    std::string reprStr;
    PyObject * dict    = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * __str__() const
  {
    PyObject * dict = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  int has_key(const std::string & name) const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * items() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * iteritems() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * iterkeys() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * itervalues() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * keys() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  void update(PyObject * dict, bool strict=true)
  {
    PyTrilinos::ResponseToIllegalParameters flag;
    if (strict) flag = PyTrilinos::raiseError;
    else        flag = PyTrilinos::storeNames;
    updateParameterListWithPyDict(dict,*self,flag);
  }

  void update(const ParameterList & plist)
  {
    self->setParameters(plist);
  }

  /******************************************************************/
  // Values method
  PyObject * values() const
  {
    PyObject * dict   = PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::ignore);
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
  PyObject * asDict() const
  {
    return PyTrilinos::parameterListToNewPyDict(*self,PyTrilinos::storeNames);
  }
}    // %extend ParameterList

%ignore Teuchos::ParameterList::set;
%ignore Teuchos::ParameterList::setEntry;
%ignore Teuchos::ParameterList::get;
%ignore Teuchos::ParameterList::getPtr;
%ignore Teuchos::ParameterList::getEntryPtr;
%ignore Teuchos::ParameterList::sublist(const std::string &) const;
%ignore Teuchos::ParameterList::isType(const std::string &) const;
%ignore Teuchos::ParameterList::isType(const std::string &, any*) const;
%ignore Teuchos::ParameterList::unused(ostream &) const;
%ignore Teuchos::ParameterList::begin() const;
%ignore Teuchos::ParameterList::end() const;
%ignore Teuchos::ParameterList::entry(ConstIterator) const;
%ignore Teuchos::ParameterList::name(ConstIterator) const;
#ifndef HAVE_TEUCHOS
%warn "HAVE_TEUCHOS IS NOT DEFINED!!!!"
#endif
%include "Teuchos_ParameterList.hpp"
// SWIG thinks that PrintOptions is an un-nested Teuchos class, so we
// need to trick the C++ compiler into understanding this so called
// un-nested Teuchos type.
%{
namespace Teuchos
{
typedef ParameterList::PrintOptions PrintOptions;
}
%}

////////////////////////////////////////////
// Teuchos::ParameterListAcceptor support //
////////////////////////////////////////////
%include "Teuchos_ParameterListAcceptor.hpp"
