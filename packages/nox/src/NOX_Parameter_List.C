// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Parameter_List.H"	// class definition

/* NOTE: ASCI Red (TFLOP) does not support the i-> funtion for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */
#if defined(TFLOP)
#define ITER (*i). 

#else
#define ITER i->

#endif

using namespace NOX::Parameter;

List::List() {}

List::List(const List& source) 
{
  params = source.params;
}

List& List::operator=(const List& source) 
{
  if (&source == this)
    return *this;

  params = source.params;
  return *this;
}

List::~List() 
{
}

void List::unused() const
{
  // Warn about any unused parameters
  /* NOTE FROM TAMMY: Note that this does not check sublists. May want
     to add that functionality later. */
  for (PCConstIterator i = params.begin(); i != params.end(); ++i) {
    if (!(ITER second.isUsed())) {
      cout << "WARNING: Parameter \"" << ITER first << "\" " << ITER second
	   << " is unused" << endl;
    }
  }
}

void List::setParameter(const string& name, bool value)
{
  params[name].setValue(value);
}

void List::setParameter(const string& name, int value)
{
  params[name].setValue(value);
}

void List::setParameter(const string& name, double value)
{
  params[name].setValue(value);
}

void List::setParameter(const string& name, const char* value)
{
  params[name].setValue(value);
}

void List::setParameter(const string& name, const string& value)
{
  params[name].setValue(value);
}


bool List::getParameter(const string& name, bool nominal)
{
  PCConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (ITER second.isBool()))
    return ITER second.getBoolValue();

  cerr << "NOX::Parameter::List::getParameter - get error for bool" << endl;
  throw "NOX Error";
}

int List::getParameter(const string& name, int nominal) 
{
  PCConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (ITER second.isInt()))
    return ITER second.getIntValue();

  cerr << "NOX::Parameter::List::getParameter - get error for int" << endl;
  throw "NOX Error";
}

double List::getParameter(const string& name, double nominal) 
{
  PCConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (ITER second.isDouble()))
    return ITER second.getDoubleValue();

  cerr << "NOX::Parameter::List::getParameter - get error for double" << endl;
  throw "NOX Error";

}

const string& List::getParameter(const string& name, const char* nominal) 
{
  PCConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (ITER second.isString()))
    return ITER second.getStringValue();

  cerr << "NOX::Parameter::List::getParameter - get error for string" << endl;
  throw "NOX Error";
}

const string& List::getParameter(const string& name, const string& nominal) 
{
  PCConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (ITER second.isString()))
    return ITER second.getStringValue();

  cerr << "NOX::Parameter::List::getParameter - get error for string" << endl;
  throw "NOX Error";
}
  
bool List::getParameter(const string& name, bool nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isBool()))
    return ITER second.getBoolValue();
  return nominal;
}

int List::getParameter(const string& name, int nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isInt()))
    return ITER second.getIntValue();
  return nominal;
}

double List::getParameter(const string& name, double nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isDouble()))
    return ITER second.getDoubleValue();
  return nominal;
}

const string& List::getParameter(const string& name, const char* nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isString()))
    return ITER second.getStringValue();

  // Save nominal char* value as a string, and return the string value.
  tmpstrings.push_back(nominal);
  return tmpstrings.back();
}

const string& List::getParameter(const string& name, const string& nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isString()))
    return ITER second.getStringValue();
  return nominal;
}
  
bool List::isParameterBool(const string& name) const
{
  PCConstIterator i = params.find(name);

  if (i != params.end())
    return (ITER second.isBool());

  return false;
}

bool List::isParameterInt(const string& name) const
{
  PCConstIterator i = params.find(name);

  if (i != params.end())
    return (ITER second.isInt());

  return false;
}

bool List::isParameterDouble(const string& name) const
{
  PCConstIterator i = params.find(name);

  if (i != params.end())
    return (ITER second.isDouble());

  return false;
}

bool List::isParameterString(const string& name) const
{
  PCConstIterator i = params.find(name);

  if (i != params.end())
    return (ITER second.isString());

  return false;
}

bool List::isParameterSublist(const string& name) const
{
  PCConstIterator i = params.find(name);

  if (i != params.end())
    return (ITER second.isList());

  return false;
}

bool List::isParameter(const string& name) const
{
  return (params.find(name) != params.end());
}

bool List::isParameterEqual(const string& name, bool value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isBool()))
    return (ITER second.getBoolValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, int value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isInt()))
    return (ITER second.getIntValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, double value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isDouble()))
    return (ITER second.getDoubleValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, const char* value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isString()))
    return (ITER second.getStringValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, const string& value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (ITER second.isString()))
    return (ITER second.getStringValue() == value);
  return false;
}


List& List::sublist(const string& name)
{
  // Find name in list, if it exists.
  PCIterator i = params.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params.end()) {
    if (ITER second.isList()) 
      return (ITER second.getListValue());
    else
      cerr << "ERROR: Parameter " << name << " is not a list." << endl;
      throw "NOX Error";
  }

  // If it does not exist, create a new empty list and return a reference
  return params[name].setList(true);
}

const List& List::sublist(const string& name) const
{
  // Find name in list, if it exists.
  PCConstIterator i = params.find(name);

  // If it does not exist, throw an error
  if (i == params.end()) {
    cerr << "ERROR: Parameter " << name << " is not a valid list." << endl;
    throw "NOX Error";
  }

  // If it does exist and is a list, return the list value.
  if (ITER second.isList()) 
    return (ITER second.getListValue());

  // Otherwise, the parameter exists but is not a list. Throw an error.
  cerr << "ERROR: Parameter " << name << " is not a list." << endl;
  throw "NOX Error";
}
  
ostream& List::print(ostream& stream, int indent) const
{
  if (params.begin() == params.end()) {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << "[empty list]" << endl;
  }
  else 
    for (PCConstIterator i = params.begin(); i != params.end(); ++i) {
      for (int j = 0; j < indent; j ++)
	stream << ' ';
      if (ITER second.isList()) {
	stream << ITER first << " -> " << endl;
	ITER second.getListValue().print(stream, indent + 2);
      }
      else
	stream << ITER first << " = " << ITER second << endl;
    }
  return stream;
}

#undef ITER
