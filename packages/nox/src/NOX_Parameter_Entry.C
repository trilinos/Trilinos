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
  isGotten(false),
  isSetByGet(false)
{
}

Entry::Entry(const Entry& source) :
  type(NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isGotten(false),
  isSetByGet(false)
{
  operator=(source);
}

Entry& Entry::operator=(const Entry& source)
{
  if (&source == this)
    return *this;

  type = source.type;
  bval = source.bval;
  ival = source.ival;
  dval = source.dval;
  sval = source.sval;
  
  delete lval;
  if ((type == LIST) && (source.lval != NULL)) {
    lval = new List(*source.lval);
  }
  
  isGotten = source.isGotten;
  isSetByGet = source.isSetByGet;
  return *this;
}

Entry::Entry(bool value, bool isCreatedByGet) : 
  type(BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet)
{
}

Entry::Entry(int value, bool isCreatedByGet) : 
  type(INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(double value, bool isCreatedByGet) : 
  type(DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const string& value, bool isCreatedByGet) : 
  type(STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::~Entry() 
{
  if (type == LIST)
    delete lval;
}

void Entry::setValue(bool value, bool isCreatedByGet)
{
  type = BOOL;
  bval = value;
  isSetByGet = isCreatedByGet;
  isGotten = false;
}

void Entry::setValue(int value, bool isCreatedByGet)
{
  type = INT;
  ival = value;
  isSetByGet = isCreatedByGet;
  isGotten = false;
}

void Entry::setValue(double value, bool isCreatedByGet)
{
  type = DOUBLE;
  dval = value;
  isSetByGet = isCreatedByGet;
  isGotten = false;
}

void Entry::setValue(const char* value, bool isCreatedByGet)
{
  type = STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
  isGotten = false;
}

void Entry::setValue(const string& value, bool isCreatedByGet)
{
  type = STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
  isGotten = false;
}

List& Entry::setList(bool isCreatedByGet)
{
  type = LIST;
  lval = new List();
  isSetByGet = isCreatedByGet;
  isGotten = true;
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
  isGotten = true;
  return bval;
}

int Entry::getIntValue() const
{
  isGotten = true;
  return ival;
}

double Entry::getDoubleValue() const
{
  isGotten = true;
  return dval;
}

const string& Entry::getStringValue() const
{
  isGotten = true;
  return sval;
}

List& Entry::getListValue() 
{
  isGotten = true;
  return *lval;
}

const List& Entry::getListValue() const
{
  isGotten = true;
  return *lval;
}

bool Entry::isUsed() const
{
  return isGotten;
}


ostream& Entry::leftshift(ostream& stream) const
{
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

  if (isSetByGet)
    cout << "   [default]";
  else if (!isGotten)
    cout << "   [unused]";
  

  return stream;
}

ostream& operator<<(ostream& stream, const Entry& e)
{
  return e.leftshift(stream);
}


