// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Parameter_Entry.H" // class definition
#include "NOX_Parameter_List.H"	 // for sublists

using namespace NOX::Parameter;

Entry::Entry() : 
  type(NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

Entry::Entry(bool value) : 
  type(BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

Entry::Entry(int value) : 
  type(INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

Entry::Entry(double value) : 
  type(DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  lval(NULL),
  isused(false) 
{
}

Entry::Entry(const string& value) : 
  type(STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  lval(NULL),
  isused(false) 
{
}

Entry::~Entry() 
{
  if (type == LIST)
    delete lval;
}

void Entry::setValue(bool value)
{
  type = BOOL;
  bval = value;
  isused = false;
}

void Entry::setValue(int value)
{
  type = INT;
  ival = value;
  isused = false;
}

void Entry::setValue(double value)
{
  type = DOUBLE;
  dval = value;
  isused = false;
}

void Entry::setValue(const char* value)
{
  type = STRING;
  sval = value;
  isused = false;
}

void Entry::setValue(const string& value)
{
  type = STRING;
  sval = value;
  isused = false;
}

List& Entry::setList()
{
  type = LIST;
  lval = new List();
  isused = true;
  return *lval;
}


bool Entry::isBool() const
{
  return (type == BOOL);
}

bool Entry::isInt() const
{
  return (type == INT);
}

bool Entry::isDouble() const
{
  return (type == DOUBLE);
}

bool Entry::isString() const
{
  return (type == STRING);
}

bool Entry::isList() const
{
  return (type == LIST);
}

bool Entry::getBoolValue() const
{
  isused = true;
  return bval;
}

int Entry::getIntValue() const
{
  isused = true;
  return ival;
}

double Entry::getDoubleValue() const
{
  isused = true;
  return dval;
}

const string& Entry::getStringValue() const
{
  isused = true;
  return sval;
}

List& Entry::getListValue() 
{
  isused = true;
  return *lval;
}

const List& Entry::getListValue() const
{
  isused = true;
  return *lval;
}

bool Entry::isUsed() const
{
  return isused;
}


ostream& Entry::leftshift(ostream& stream) const
{
  string sep = "     ";
  switch(type) {
  case BOOL: 
    stream << (bval ? "true" : "false");
    break;
  case INT:
    stream << ival;
    break;
  case DOUBLE:
    stream << dval;
    break;
  case STRING:
    stream << "\"" << sval << "\"";
    break;
  case LIST:
    break;
  default:
    stream << "(empty non-typed parameter)";
    break;
  }
  return stream;
}

ostream& operator<<(ostream& stream, const Entry& e)
{
  return e.leftshift(stream);
}


