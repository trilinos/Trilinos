// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_ParameterList.H"	// class definition

//----------------------------------------------------------------------
// NLS_ParameterList
//----------------------------------------------------------------------

NLS_ParameterList::NLS_ParameterList() {}

NLS_ParameterList::~NLS_ParameterList() 
{
}

void NLS_ParameterList::unused() const
{
  // Warn about any unused parameters
  /* NOTE FROM TAMMY: Note that this does not check sublists. May want
     to add that functionality later. */
  for (PCConstIterator i = params.begin(); i != params.end(); ++i) {
    if (!(i->second.isUsed())) {
      cout << "WARNING: Parameter \"" << i->first << "\" " << i->second
	   << " is unused" << endl;
    }
  }
}

void NLS_ParameterList::setParameter(const string& name, bool value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, int value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, double value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, const char* value)
{
  params[name].setValue(value);
}

void NLS_ParameterList::setParameter(const string& name, const string& value)
{
  params[name].setValue(value);
}


bool NLS_ParameterList::getParameter(const string& name, bool nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isBool()))
    return i->second.getBoolValue();
  return nominal;
}

int NLS_ParameterList::getParameter(const string& name, int nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isInt()))
    return i->second.getIntValue();
  return nominal;
}

double NLS_ParameterList::getParameter(const string& name, double nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isDouble()))
    return i->second.getDoubleValue();
  return nominal;
}

const string& NLS_ParameterList::getParameter(const string& name, const char* nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return i->second.getStringValue();

  // Save nominal char* value as a string, and return the string value.
  tmpstrings.push_back(nominal);
  return tmpstrings[tmpstrings.size() - 1];
}

const string& NLS_ParameterList::getParameter(const string& name, const string& nominal) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return i->second.getStringValue();
  return nominal;
}
  
bool NLS_ParameterList::isParameter(const string& name) const
{
  return (params.find(name) != params.end());
}

bool NLS_ParameterList::isParameterEqual(const string& name, bool value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isBool()))
    return (i->second.getBoolValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, int value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isInt()))
    return (i->second.getIntValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, double value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isDouble()))
    return (i->second.getDoubleValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, const char* value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return (i->second.getStringValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, const string& value) const
{
  PCConstIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return (i->second.getStringValue() == value);
  return false;
}


NLS_ParameterList& NLS_ParameterList::sublist(const string& name)
{
  // Find name in list, if it exists.
  PCIterator i = params.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params.end()) {
    if (i->second.isList()) 
      return (i->second.getListValue());
    else
      cerr << "ERROR: Parameter " << name << " is not a list." << endl;
      throw;
  }

  // If it does not exist, create a new empty list and return a reference
  return params[name].setList();
}
  
ostream& NLS_ParameterList::print(ostream& stream, int indent = 0) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "-- List Entries --\n";
  for (PCConstIterator i = params.begin(); i != params.end(); ++i) {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << i->first << " = " << i->second << endl;
    if (i->second.isList()) 
      i->second.getListValue().print(stream, indent + 2);
  }
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "-- End Of List --\n";
  return stream;
}
