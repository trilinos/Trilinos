// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _VARIABLERTC_H
#define _VARIABLERTC_H

#include "RTC_commonRTC.hh"
#include "RTC_ValueRTC.hh"

#include <cassert>
#include <string>

namespace PG_RuntimeCompiler {

/**
 * A Varaible object represents the variables in the code the user gives us.
 */

class Variable: public Value
{
 public:

  /**
   * Constructor -> Trivial
   *
   * @param name    - The name of the variable
   * @param type    - The data type of the variable
   * @param objType - The object type of the variable
   */
  Variable(const std::string& name, Type type, ObjectType objType)
    : Value(type, objType) {
    _name    = name;
    _address = NULL;

    _willBeInitAtTimeOfUse = false;
  }

  /**
   * getName -> This method returns the name of the variable
   */
  std::string getName() const {return _name;}

  /**
   * setAddress -> This method sets the memory address of a variable. This
   *               only applies to Variables that have been passed in to
   *               Function by-reference.
   *
   * @param addr - The new address
   */
  virtual void setAddress(void* addr) = 0;

  /**
   * setSize - Applies to arrays only, no-op here
   */
  virtual void setSize(long /* size */) {}

  /**
   * evaluateSizeExpr - Applies to arrays only, no-op here
   */
  virtual void evaluateSizeExpr() {}

  /**
   * init -> This method sets _willBeInitAtTimeOfUse to true
   */
  void init() { _willBeInitAtTimeOfUse = true;}

  /**
   * isInit -> This method returns _willBeInitAtTimeOfUse
   */
  bool isInit() const { return _willBeInitAtTimeOfUse;}

 protected:

  std::string _name; //!< The name of the variable

  bool _willBeInitAtTimeOfUse; /**!< helps us find errors where the user is
                                *    trying to use an uninitialized variable.
                                */

  void* _address; //!< The address location of the variable
};

}

#endif
