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
  isused(false) 
{
}

Entry::Entry(const Entry& source) :
  type(NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  lval(NULL),
  isused(false) 
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
  
  isused = source.isused;
  return *this;
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


