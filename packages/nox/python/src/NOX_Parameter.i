// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package="PyTrilinos.NOX") Parameter

%{
// System includes
#include <sstream>
#include <Python.h>

// NOX includes
#include "NOX_Parameter_List.H"
#include "Utils_enums.H"
%}

// Ignore directives
%ignore NOX::Parameter::List::operator=(const List&);
%ignore *::print(ostream& stream, int indent = 0) const;
%ignore NOX::Parameter::List::getParameter(string const &, bool        	 );
%ignore NOX::Parameter::List::getParameter(string const &, int         	 );
%ignore NOX::Parameter::List::getParameter(string const &, double      	 );
%ignore NOX::Parameter::List::getParameter(string const &, string      	 );
%ignore NOX::Parameter::List::getParameter(string const &, const char *	 );
%ignore NOX::Parameter::List::getParameter(string const &, bool        	 ) const;
%ignore NOX::Parameter::List::getParameter(string const &, int         	 ) const;
%ignore NOX::Parameter::List::getParameter(string const &, double      	 ) const;
%ignore NOX::Parameter::List::getParameter(string const &, string      	 ) const;
%ignore NOX::Parameter::List::getParameter(string const &, string const &) const;
%ignore NOX::Parameter::List::getParameter(string const &, const char *  ) const;
%ignore NOX::Parameter::List::getParameter(string const &, Arbitrary const &) const;
%ignore NOX::Parameter::List::isParameterEqual(string const &, const char *) const;
%ignore NOX::Parameter::List::setParameter(string const &, const char *);
%ignore NOX::Parameter::List::setParameter(string const &, char *      );
%ignore NOX::Parameter::List::sublist(     string const &              ) const;

// Rename directives
/* None */

// Define new int and double typecheck typemaps to support bools
// The precedence of bools was changed to be after int's and double's
// and these ensure that bools don't match the int and double typechecks
namespace NOX {
  namespace Parameter {
    %typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) int {
      $1 = PyInt_CheckExact($input) ? 1 : 0;
    }
    %typemap(typecheck,precedence=SWIG_TYPECHECK_FLOAT) double {
      $1 = PyFloat_CheckExact($input) ? 1 : 0;
    }
  }
}

// SWIG library includes
%include "std_string.i"

// NOX interface imports
%import "Utils_enums.H"

// NOX interface includes
using namespace std;
%include "NOX_Parameter_List.H"

// Extensions
%extend NOX::Parameter::List {

  PyObject * getParameter(const string & name) {
    if (self->isParameterInt(name)) {
      int value(0);
      return PyInt_FromLong(long(self->getParameter(name, value)));
    } else if (self->isParameterBool(name)) {
      bool value(0);
      return PyBool_FromLong(self->getParameter(name, value));
    } else if (self->isParameterDouble(name)) {
      double value(0.0);
      return PyFloat_FromDouble(self->getParameter(name, value));
    } else if (self->isParameterString(name)) {
      string value("");
      return PyString_FromString(self->getParameter(name, value).c_str());
    } else {
      PyErr_SetString(PyExc_TypeError,"Unknown Parameter Type");
      return NULL;
    }
  }

  string __str__() {
     stringstream os;
     self->print(os);                  // Put the output in os
     string s = os.str();              // Extract the string from os
     return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}
