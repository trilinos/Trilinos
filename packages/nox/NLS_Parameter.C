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
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(bool value) : 
  type(NONE),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(int value) : 
  type(NONE),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(double value) : 
  type(NONE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  isused(false) 
{
}

NLS_Parameter::NLS_Parameter(const string& value) : 
  type(NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  isused(false) 
{
}

NLS_Parameter::~NLS_Parameter() 
{
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

void NLS_Parameter::setValue(const string& value)
{
  type = STRING;
  sval = value;
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


bool NLS_Parameter::getBoolValue()
{
  return bval;
}

int NLS_Parameter::getIntValue()
{
  return ival;
}

double NLS_Parameter::getDoubleValue()
{
  return dval;
}

string& NLS_Parameter::getStringValue()
{
  return sval;
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
    stream << "boolean" << sep << bval;
    break;
  case INT:
    stream << "integer" << sep << ival;
    break;
  case DOUBLE:
    stream << "double" << sep << dval;
    break;
  case STRING:
    stream << "string" << sep << sval;
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


