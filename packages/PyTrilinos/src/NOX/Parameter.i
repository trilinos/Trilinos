// -*- c++ -*-

%module Parameter

%{
// System includes
#include <sstream>
#include <Python.h>

// NOX includes
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

// Namespace flattening
using namespace NOX           ;
using namespace NOX::Parameter;
%}

// Ignore directives
%ignore NOX::Parameter::List::operator=(const List&);
%ignore *::print(ostream &, int) const;
%ignore NOX::Parameter::List::getParameter(string const &, bool        );
%ignore NOX::Parameter::List::getParameter(string const &, int         );
%ignore NOX::Parameter::List::getParameter(string const &, double      );
%ignore NOX::Parameter::List::getParameter(string const &, string      );
%ignore NOX::Parameter::List::getParameter(string const &, const char *);
%ignore NOX::Parameter::List::getParameter(string const &, bool        ) const;
%ignore NOX::Parameter::List::getParameter(string const &, int         ) const;
%ignore NOX::Parameter::List::getParameter(string const &, double      ) const;
%ignore NOX::Parameter::List::getParameter(string const &, string      ) const;
%ignore NOX::Parameter::List::getParameter(string const &, const char *) const;
%ignore NOX::parameter::List::isParameterBool(const string &) const;
%ignore NOX::Parameter::List::isParameterEqual(string const &, bool        ) const;
%ignore NOX::Parameter::List::isParameterEqual(string const &, const char *) const;
%ignore NOX::Parameter::List::setParameter(string const &, bool        );
%ignore NOX::Parameter::List::setParameter(string const &, const char *);
%ignore NOX::Parameter::List::setParameter(string const &, char *      );

// Rename directives
/* None */

// SWIG library includes
%include "std_string.i"

// NOX interface includes
using namespace std;
%include "NOX_Parameter_List.H"

// Extensions
%extend NOX::Parameter::List {

  PyObject * getParameter(const string & name) {
    if (self->isParameterInt(name)) {
      int value;
      return PyInt_FromLong(long(self->getParameter(name, value)));
    } else if (self->isParameterDouble(name)) {
      double value;
      return PyFloat_FromDouble(self->getParameter(name, value));
    } else if (self->isParameterString(name)) {
      string value;
      return PyString_FromString(self->getParameter(name, value).c_str());
    } else {
      PyErr_SetString(PyExc_TypeError,"Unknown Parameter Type");
      return NULL;
    }
  }

  const PyObject * getParameter(const string & name) const {
    if (self->isParameterInt(name)) {
      int value;
      return PyInt_FromLong(long(self->getParameter(name, value)));
    } else if (self->isParameterDouble(name)) {
      double value;
      return PyFloat_FromDouble(self->getParameter(name, value));
    } else if (self->isParameterString(name)) {
      string value;
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

class Utils {
public:
  enum MsgType { 
    Error                    = NOX::Utils::Error,
    Warning                  = NOX::Utils::Warning,
    OuterIteration           = NOX::Utils::OuterIteration,
    InnerIteration           = NOX::Utils::InnerIteration,
    Parameters               = NOX::Utils::Parameters,
    Details                  = NOX::Utils::Details,
    OuterIterationStatusTest = NOX::Utils::OuterIterationStatusTest,
    LinearSolverDetails      = NOX::Utils::LinearSolverDetails
  };
};
