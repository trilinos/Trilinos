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

#include "Teuchos_ParameterList.hpp" 

using namespace AMESOS::Parameter;

Entry::Entry() : 
  type(AMESOS_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  aval(NULL),
  lval(NULL),
  isGotten(false),
  isSetByGet(false)
{
}

Entry::Entry(const Entry& source) :
  type(AMESOS_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  aval(NULL),
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

  reset();

  type = source.type;
  bval = source.bval;
  ival = source.ival;
  dval = source.dval;
  sval = source.sval;
  
  if ((type == AMESOS_ARBITRARY) && (source.aval != NULL)) {
    aval = source.aval->clone();
  }
  
  if ((type == AMESOS_LIST) && (source.lval != NULL)) {
    lval = new List(*source.lval);
  }
  
  isGotten = source.isGotten;
  isSetByGet = source.isSetByGet;

  return *this;
}

Entry::Entry(bool value, bool isCreatedByGet) : 
  type(AMESOS_BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""), 
  aval(NULL),
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet)
{
}

Entry::Entry(int value, bool isCreatedByGet) : 
  type(AMESOS_INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""), 
  aval(NULL),
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(double value, bool isCreatedByGet) : 
  type(AMESOS_DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""), 
  aval(NULL),
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const string& value, bool isCreatedByGet) : 
  type(AMESOS_STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value), 
  aval(NULL),
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::Entry(const Arbitrary& value, bool isCreatedByGet) : 
  type(AMESOS_STRING),
  bval(false),
  ival(0),
  dval(0),
  sval("" ),
  aval(value.clone()),
  lval(NULL),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

Entry::~Entry() 
{
  reset();
}

void Entry::reset()
{
  type = AMESOS_NONE;

  delete aval;
  aval = NULL;

  delete lval;
  lval = NULL;

  isGotten = false;
  isSetByGet = false;
}

void Entry::setValue(bool value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_BOOL;
  bval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(int value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_INT;
  ival = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(double value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_DOUBLE;
  dval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const char* value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const string& value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void Entry::setValue(const Arbitrary& value, bool isCreatedByGet)
{
  reset();
  type = AMESOS_ARBITRARY;
  aval = value.clone();
  isSetByGet = isCreatedByGet;
}

List& Entry::setList(bool isCreatedByGet)
{
  reset();
  type = AMESOS_LIST;
  lval = new List();
  isSetByGet = isCreatedByGet;
  isGotten = true;
  return *lval;
}


bool Entry::isBool() const
{
  return (type == AMESOS_BOOL);
}

bool Entry::isInt() const
{
  return (type == AMESOS_INT);
}

bool Entry::isDouble() const
{
  return (type == AMESOS_DOUBLE);
}

bool Entry::isString() const
{
  return (type == AMESOS_STRING);
}

bool Entry::isArbitrary() const
{
  return (type == AMESOS_ARBITRARY);
}

bool Entry::isList() const
{
  return (type == AMESOS_LIST);
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

const Arbitrary& Entry::getArbitraryValue() const
{
  isGotten = true;
  return *aval;
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
  case AMESOS_BOOL: 
    stream << (bval ? "true" : "false");
    break;
  case AMESOS_INT:
    stream << ival;
    break;
  case AMESOS_DOUBLE:
    stream << dval;
    break;
  case AMESOS_STRING:
    stream << "\"" << sval << "\"";
    break;
  case AMESOS_ARBITRARY:
    stream << aval->getType();
    break;
  case AMESOS_LIST:
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


#endif
