// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_ParameterList.H"

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
  for (PCIterator i = params.begin(); i != params.end(); ++i) {
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

//PRIVATE
bool NLS_ParameterList::isRecursive(const NLS_ParameterList& l) const
{
  // Check if l or any of its sublists points to "this"
  if (&l == this)
    return true;

  for (PCIterator i = l.params.begin(); i != l.params.end(); ++i) {
    if ((i->second.isList()) && 
	(isRecursive(i->second.getListValue())))
      return true;
  }

  return false;
}

bool NLS_ParameterList::setParameter(const string& name, const NLS_ParameterList& value)
{
  // Cannot add this list if it or any of its sublists is this!
  if (isRecursive(value)) 
    return false;

  params[name].setValue(value);
  return true;
}


bool NLS_ParameterList::getParameter(const string& name, bool nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isBool()))
    return i->second.getBoolValue();
  return nominal;
}

int NLS_ParameterList::getParameter(const string& name, int nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isInt()))
    return i->second.getIntValue();
  return nominal;
}

double NLS_ParameterList::getParameter(const string& name, double nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isDouble()))
    return i->second.getDoubleValue();
  return nominal;
}

const string& NLS_ParameterList::getParameter(const string& name, const char* nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return i->second.getStringValue();

  // Save nominal char* value as a string, and return the string value.
  tmpstrings.push_back(nominal);
  return tmpstrings[tmpstrings.size() - 1];
}

const string& NLS_ParameterList::getParameter(const string& name, const string& nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return i->second.getStringValue();
  return nominal;
}
  
const NLS_ParameterList& NLS_ParameterList::getParameter(const string& name, const NLS_ParameterList& nominal) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isList()))
    return i->second.getListValue();
  return nominal;
}
  

bool NLS_ParameterList::isParameter(const string& name) const
{
  return (params.find(name) != params.end());
}

bool NLS_ParameterList::isParameterEqual(const string& name, bool value) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isBool()))
    return (i->second.getBoolValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, int value) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isInt()))
    return (i->second.getIntValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, double value) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isDouble()))
    return (i->second.getDoubleValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, const char* value) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return (i->second.getStringValue() == value);
  return false;
}

bool NLS_ParameterList::isParameterEqual(const string& name, const string& value) const
{
  PCIterator i = params.find(name);
  if ((i != params.end()) && (i->second.isString()))
    return (i->second.getStringValue() == value);
  return false;
}


ostream& NLS_ParameterList::print(ostream& stream, int indent = 0) const
{
  for (PCIterator i = params.begin(); i != params.end(); ++i) {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << i->first << " = " << i->second << endl;
    if (i->second.isList()) 
      i->second.getListValue().print(stream, indent + 2);
  }
  return stream;
}
