// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Parameter.H"

NLS_Parameter::NLS_Parameter() : 
  type(NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(bool value) : 
  type(BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(int value) : 
  type(INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(double value) : 
  type(DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(const string& value) : 
  type(STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  lval(NULL),
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(const NLS_ParameterList& value) : 
  type(LIST),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(&value),
  isused(false) 
{
}

NLS_Parameter::~NLS_Parameter() 
{
  // Do not delete list - the calling program is responsible for that.
}

void NLS_Parameter::setValue(bool value)
{
  type = BOOL;
  bval = value;
  isused = false;
}

void NLS_Parameter::setValue(int value)
{
  type = INT;
  ival = value;
  isused = false;
}

void NLS_Parameter::setValue(double value)
{
  type = DOUBLE;
  dval = value;
  isused = false;
}

void NLS_Parameter::setValue(const char* value)
{
  type = STRING;
  sval = value;
  isused = false;
}

void NLS_Parameter::setValue(const string& value)
{
  type = STRING;
  sval = value;
  isused = false;
}

void NLS_Parameter::setValue(const NLS_ParameterList& value)
{
  type = LIST;
  lval = &value;
  isused = false;
}


bool NLS_Parameter::isBool() const
{
  return (type == BOOL);
}

bool NLS_Parameter::isInt() const
{
  return (type == INT);
}

bool NLS_Parameter::isDouble() const
{
  return (type == DOUBLE);
}

bool NLS_Parameter::isString() const
{
  return (type == STRING);
}

bool NLS_Parameter::isList() const
{
  return (type == LIST);
}

bool NLS_Parameter::getBoolValue() const
{
  isused = true;
  return bval;
}

int NLS_Parameter::getIntValue() const
{
  isused = true;
  return ival;
}

double NLS_Parameter::getDoubleValue() const
{
  isused = true;
  return dval;
}

const string& NLS_Parameter::getStringValue() const
{
  isused = true;
  return sval;
}

const NLS_ParameterList& NLS_Parameter::getListValue() const
{
  isused = true;
  return *lval;
}

bool NLS_Parameter::isUsed() const
{
  return isused;
}


ostream& NLS_Parameter::leftshift(ostream& stream) const
{
  const string sep = ",";
  stream << "<"; 
  switch(type) {
  case BOOL: 
    stream << "boolean" << sep << (bval ? "true" : "false");
    break;
  case INT:
    stream << "integer" << sep << ival;
    break;
  case DOUBLE:
    stream << "double" << sep << dval;
    break;
  case STRING:
    stream << "string" << sep << "\"" << sval << "\"";
    break;
  case LIST:
    stream << "sublist";
    break;
  default:
    stream << "NONE";
    break;
  }
  stream << ">";
  return stream;
}

ostream& operator<<(ostream& stream, const NLS_Parameter& e)
{
  return e.leftshift(stream);
}


