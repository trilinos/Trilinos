// $Id$ 
// $Source$ 

// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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

#ifdef OLDLIST

#include "Amesos_Parameter_List.h"	// class definition

/* NOTE: ASCI Red (TFLOP) does not support the i-> funtion for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */
using namespace AMESOS::Parameter;



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
  for (ConstIterator i = params.begin(); i != params.end(); ++i) {
    if (!(entry(i).isUsed())) {
      cout << "WARNING: Parameter \"" << name(i) << "\" " << entry(i)
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

void List::setParameter(const string& name, const Arbitrary& value)
{
  params[name].setValue(value);
}


bool List::getParameter(const string& name, bool nominal)
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isBool()))
    return entry(i).getBoolValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for bool" << endl;
  throw "AMESOS Error";
}

int List::getParameter(const string& name, int nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isInt()))
    return entry(i).getIntValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for int" << endl;
  throw "AMESOS Error";
}

double List::getParameter(const string& name, double nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isDouble()))
    return entry(i).getDoubleValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for double" << endl;
  throw "AMESOS Error";

}

const string& List::getParameter(const string& name, const char* nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for string" << endl;
  throw "AMESOS Error";
}

const string& List::getParameter(const string& name, const string& nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for string" << endl;
  throw "AMESOS Error";
}
  
const Arbitrary& List::getParameter(const string& name, const Arbitrary& nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isArbitrary()))
    return entry(i).getArbitraryValue();

  cerr << "AMESOS::Parameter::List::getParameter - get error for arbitrary parameter" << endl;
  throw "AMESOS Error";
}
  
bool List::getParameter(const string& name, bool nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isBool()))
    return entry(i).getBoolValue();
  return nominal;
}

int List::getParameter(const string& name, int nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isInt()))
    return entry(i).getIntValue();
  return nominal;
}

double List::getParameter(const string& name, double nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isDouble()))
    return entry(i).getDoubleValue();
  return nominal;
}

const string& List::getParameter(const string& name, const char* nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  // Save nominal char* value as a string, and return the string value.
  tmpstrings.push_back(nominal);
  return tmpstrings.back();
}

const string& List::getParameter(const string& name, const string& nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();
  return nominal;
}
  
const Arbitrary& List::getParameter(const string& name, const Arbitrary& nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isArbitrary()))
    return entry(i).getArbitraryValue();
  return nominal;
}
  
const Arbitrary& List::getArbitraryParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isArbitrary()))
    return entry(i).getArbitraryValue();

  cerr << "AMESOS::Parameter::List::getArbitraryParameter - no such parameter" << endl;
  throw "AMESOS Error";
}
  
bool List::isParameterBool(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isBool());

  return false;
}

bool List::isParameterInt(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isInt());

  return false;
}

bool List::isParameterDouble(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isDouble());

  return false;
}

bool List::isParameterString(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isString());

  return false;
}

bool List::isParameterArbitrary(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isArbitrary());

  return false;
}

bool List::isParameterSublist(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isList());

  return false;
}

bool List::isParameter(const string& name) const
{
  return (params.find(name) != params.end());
}

bool List::isParameterEqual(const string& name, bool value) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isBool()))
    return (entry(i).getBoolValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, int value) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isInt()))
    return (entry(i).getIntValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, double value) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isDouble()))
    return (entry(i).getDoubleValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, const char* value) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return (entry(i).getStringValue() == value);
  return false;
}

bool List::isParameterEqual(const string& name, const string& value) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return (entry(i).getStringValue() == value);
  return false;
}


List& List::sublist(const string& name)
{
  // Find name in list, if it exists.
  Iterator i = params.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params.end()) {
    if (entry(i).isList()) 
      return (entry(i).getListValue());
    else
      cerr << "ERROR: Parameter " << name << " is not a list." << endl;
      throw "AMESOS Error";
  }

  // If it does not exist, create a new empty list and return a reference
  return params[name].setList(true);
}

const List& List::sublist(const string& name) const
{
  // Find name in list, if it exists.
  ConstIterator i = params.find(name);

  // If it does not exist, throw an error
  if (i == params.end()) {
    cerr << "ERROR: Parameter " << name << " is not a valid list." << endl;
    throw "AMESOS Error";
  }

  // If it does exist and is a list, return the list value.
  if (entry(i).isList()) 
    return (entry(i).getListValue());

  // Otherwise, the parameter exists but is not a list. Throw an error.
  cerr << "ERROR: Parameter " << name << " is not a list." << endl;
  throw "AMESOS Error";
}
  
ostream& List::print(ostream& stream, int indent) const
{
  if (params.begin() == params.end()) 
  {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << "[empty list]" << endl;
  }
  else 
    for (ConstIterator i = params.begin(); i != params.end(); ++i) 
    {
      for (int j = 0; j < indent; j ++)
	stream << ' ';
      if (entry(i).isList()) 
      {
	stream << name(i) << " -> " << endl;
	entry(i).getListValue().print(stream, indent + 2);
      }
      else if (entry(i).isArbitrary()) 
      {
	stream << name(i) << " = " << entry(i) << endl;
	entry(i).getArbitraryValue().print(stream, indent + 2);
      }
      else
	stream << name(i) << " = " << entry(i) << endl;
    }
  return stream;
}


#if defined(TFLOP)

const string& List::name(ConstIterator i) const
{
  return ((*i).first);
}

Entry& List::entry(Iterator i)
{
  return ((*i).second);
}

const Entry& List::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else

const string& List::name(ConstIterator i) const
{
  return (i->first);
}

Entry& List::entry(Iterator i)
{
  return (i->second);
}

const Entry& List::entry(ConstIterator i) const
{
  return (i->second);
}

#endif

#endif
