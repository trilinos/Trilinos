// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _VALUERTC_H
#define _VALUERTC_H

#include "RTC_ObjectRTC.hh"
#include "RTC_commonRTC.hh"

#include <string>
#include <cassert>

namespace PG_RuntimeCompiler {

/**
 * A Value object represents the operands in the code the user gives us.
 */

class Value: public Object
{
 protected:
  Type _type; //The data type of the Value

 public:

  /**
   * Constructor -> Trivial
   *
   * @param type    - The data type of the value
   * @param objType - The object type of the value
   */
  Value(Type type, ObjectType objType) : Object(objType) {_type = type;}

  /**
   * getType -> This method returns the Type of the Value.
   */
  Type getType() const {return _type;}

  /**
   * getValue -> Should not be called
   */
  virtual double getValue() { assert(false); return 0;}

  /**
   * getValue -> Should not be called
   */
  virtual double getArrayValue(long /* offset */) const { assert(false); return 0;}

  /**
   * setValue -> Should not be called
   */
  virtual void setValue(double /* value */) { assert(false); }

  /**
   * setArrayValue -> Should not be called
   */
  virtual void setArrayValue(double /* value */, long /* offset */) { assert(false); }

  /**
   * getSize - The size of everything is zero, except for arrays
   */
  virtual long getSize() const { return 0;}
};

}

#endif
