// -*- c++ -*-

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

%{
// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"
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
%#if PY_VERSION_HEX >= 0x03000000
    PyObject * byteObj = PyUnicode_AsASCIIString(strObj);
    if (!byteObj) throw PyTrilinos::PythonException();
    self->addAttribute(name, std::string(PyBytes_AsString(byteObj)));
    Py_DECREF(byteObj);
%#else
    self->addAttribute(name, std::string(PyString_AsString(strObj)));
%#endif
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
