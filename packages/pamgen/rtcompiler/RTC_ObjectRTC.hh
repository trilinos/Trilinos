// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _OBJECTRTC_H
#define _OBJECTRTC_H

#include "RTC_commonRTC.hh"

#include <iostream>

namespace PG_RuntimeCompiler {

/**
 * A Object object represnts the operands and operators in the code the user
 * gives us. It can be a Variable, constant (Number), ArrayIndex, etc. The
 * only thing all object have in common is they know what type of object they
 * are.
 */

class Object
{
 public:

  /**
   * Constructor -> Sets the object's type
   *
   * @param type - The type of the object being created
   */
  Object(ObjectType type) { _objType = type;}

  /**
   * Destructor -> The destructor is a no-op
   */
  virtual ~Object() {}

  virtual std::ostream& operator<<(std::ostream& os) const = 0;

  /**
   * getObjectType -> Returns the type of the object
   */
  ObjectType getObjectType() const { return _objType;}

 protected:
  ObjectType _objType; //!< The type of the object
};

std::ostream& operator<<(std::ostream& os, const Object& obj);

}

#endif
