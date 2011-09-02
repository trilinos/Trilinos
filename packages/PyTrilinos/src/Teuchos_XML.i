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
#include "Teuchos_XMLObjectImplem.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StringInputSource.hpp"

// PyTrilinos includes
#include "PyTrilinos_PythonException.h"
%}

// Teuchos imports
%import "Teuchos_XMLObjectImplem.hpp"

////////////////////////////////
// Teuchos::XMLObject support //
////////////////////////////////
%feature("docstring")
Teuchos::XMLObject
"Representation of an XML data tree."
%feature("docstring")
Teuchos::XMLObject::XMLObject
"The constructor that takes an ``XMLObjectImplem*`` argument has been
removed.  The ``XMLObjectImplem`` class is hidden from the python
user."
%feature("docstring")
Teuchos::XMLObject::addAttribute
"The ``addAttribute(name, value)`` method converts the value to its
string representation and adds it as an attribute to the XML object."
%feature("docstring")
Teuchos::XMLObject::getWithDefault
"The ``getWithDefault(name, defaultValue)`` method will return the
typed value of attribute ``name``.  If the attribute does not exist,
then ``defaultValue``, which can be any python object, is returned.
If the underlying object cannot be converted to a python object, then
a python string is returned."
%feature("docstring")
Teuchos::XMLObject::__str__
"The ``__str__()`` method is provided so that it is possible to
``print`` an ``XMLObject`` object.  It returns the same string as the
``toString()`` method, but if ``toString()`` raises an exception (such
as when the ``XMLObject`` is empty), the ``__str__()`` method returns
the empty string."
%ignore Teuchos::XMLObject::XMLObject(XMLObjectImplem*);
%extend Teuchos::XMLObject
{
  void addAttribute(const std::string & name, PyObject * value)
  {
    PyObject * strObj = PyObject_Str(value);
    if (!strObj) throw PyTrilinos::PythonException();
    self->addAttribute(name, std::string(PyString_AsString(strObj)));
    Py_DECREF(strObj);
  }
  PyObject * getWithDefault(const std::string & name, PyObject * defaultValue)
  {
    PyObject * value = defaultValue;
    if (self->hasAttribute(name))
    {
      PyObject * globals = PyDict_New();
      const char * attrStr = self->getAttribute(name).c_str();
      value = PyRun_String(attrStr, Py_file_input, globals, globals);
      if (!value)
      {
	PyErr_Clear();
	value = PyString_FromString(attrStr);
      }
      Py_DECREF(globals);
      if (!value) throw PyTrilinos::PythonException();
    }
    return value;
  }
  std::string __str__() const
  {
    try
    {
      return self->toString();
    }
    catch(std::logic_error e)
    {
      return std::string("");
    }
  }
}
%include "Teuchos_XMLObject.hpp"

/////////////////////////////////////////////
// Teuchos::XMLParameterListReader support //
/////////////////////////////////////////////
namespace Teuchos
{
%ignore XMLParameterListReader::toParameterList(const XMLObject&,
						RCP< DependencySheet >) const;
}
%include "Teuchos_XMLParameterListReader.hpp"

/////////////////////////////////////////////
// Teuchos::XMLParameterListWriter support //
/////////////////////////////////////////////
namespace Teuchos
{
%extend XMLParameterListWriter
{
  XMLObject toXML(const ParameterList & p) const
  {
    return self->toXML(p);
  }
}
%ignore XMLParameterListWriter::toXML(const ParameterList&,
				      RCP< const DependencySheet >) const;
}
%include "Teuchos_XMLParameterListWriter.hpp"

/////////////////////////////////////
// Teuchos::XMLInputSource support //
/////////////////////////////////////
%ignore Teuchos::XMLInputSource::stream() const;
%include "Teuchos_XMLInputSource.hpp"

//////////////////////////////////////
// Teuchos::FileInputSource support //
//////////////////////////////////////
%ignore Teuchos::FileInputSource::stream() const;
%include "Teuchos_FileInputSource.hpp"

////////////////////////////////////////
// Teuchos::StringInputSource support //
////////////////////////////////////////
%ignore Teuchos::StringInputSource::stream() const;
%include "Teuchos_StringInputSource.hpp"
