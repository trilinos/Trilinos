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

// A useful typedef for manipulating the params map
typedef map<string, NLS_Parameter>::const_iterator paramsiter;

NLS_ParameterList::NLS_ParameterList() {}

NLS_ParameterList::~NLS_ParameterList() 
{
}

void NLS_ParameterList::unused() const
{
  // Warn about any unused parameters
  for (paramsiter it = params.begin(); it != params.end(); it ++) {
    if (!(it->second.isUsed())) {
      cout << "WARNING: Parameter \"" << it->first << "\" " << it->second
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

bool NLS_ParameterList::isRecursive(const List& l) const
{
  // Check if l or any of its sublists points to "this"
  if (&l == this)
    return true;

  for (paramsiter it = l.params.begin(); it != l.params.end(); it ++) {
    if ((it->second.isList()) && 
	(isRecursive(it->second.getListValue())))
      return true;
  }

  return false;
}

bool NLS_ParameterList::setParameter(const string& name, const List& value)
{
  // Cannot add this list if it or any of its sublists is this!
  if (isRecursive(value)) 
    return false;

  params[name].setValue(value);
  return true;
}


bool NLS_ParameterList::getParameter(const string& name, bool nominal) const
{
  paramsiter it = params.find(name);
  if ((it != params.end()) && (it->second.isBool()))
    return it->second.getBoolValue();
  return nominal;
}

int NLS_ParameterList::getParameter(const string& name, int nominal) const
{
  paramsiter it = params.find(name);
  if ((it != params.end()) && (it->second.isInt()))
    return it->second.getIntValue();
  return nominal;
}

double NLS_ParameterList::getParameter(const string& name, double nominal) const
{
  paramsiter it = params.find(name);
  if ((it != params.end()) && (it->second.isDouble()))
    return it->second.getDoubleValue();
  return nominal;
}

const string& NLS_ParameterList::getParameter(const string& name, const char* nominal) const
{
  tmpstrings.push_back(nominal);
  return getParameter(name, tmpstrings[tmpstrings.size() - 1]);
}

const string& NLS_ParameterList::getParameter(const string& name, const string& nominal) const
{
  paramsiter it = params.find(name);
  if ((it != params.end()) && (it->second.isString()))
    return it->second.getStringValue();
  return nominal;
}
  
const List& NLS_ParameterList::getParameter(const string& name, const List& nominal) const
{
  paramsiter it = params.find(name);
  if ((it != params.end()) && (it->second.isList()))
    return it->second.getListValue();
  return nominal;
}
  

ostream& NLS_ParameterList::print(ostream& stream, int indent = 0) const
{
  for (paramsiter it = params.begin(); it != params.end(); it ++) {
    for (int i = 0; i < indent; i ++)
      stream << ' ';
    stream << it->first << " = " << it->second << endl;
    if (it->second.isList()) 
      it->second.getListValue().print(stream, indent + 2);
  }
  return stream;
}
