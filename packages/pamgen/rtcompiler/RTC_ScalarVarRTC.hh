// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _SCALARVAR_HH
#define _SCALARVAR_HH

#include "RTC_commonRTC.hh"
#include "RTC_ValueRTC.hh"

#include <cassert>
#include <string>
#include <iostream>

namespace PG_RuntimeCompiler {

/**
 * ScalarVar object will represent plain, non-array variables
 */

template <class T>
class ScalarVar : public Variable
{
 private:
  T _value; //!< The value contained by the variable

 public:

  /**
   * Constructor -> Trivial
   *
   * @param name  - The name of the variable
   * @param value - The initial value of the variable (optional)
   */
  ScalarVar(const std::string& name, T value = 0)
    : Variable(name, TypeToTypeT<T>::value, ScalarVarOT)
  {
    _value = value;
  }

  std::ostream& operator<<(std::ostream& os) const {
    os << _value;
    return os;
  }

  /**
   * setValue -> Changes the value of a variable.
   *
   * @param value - The variable's new value
   */
  void setValue(double value) {
    _value = (T) value;

    if (_address != NULL)
      *((T*)_address) = _value;
  }

  /**
   * setAddress -> Sets the addres of the variable. Used by varAddrFill method.
   *
   * @param addr - The variable's new address
   */
  void setAddress(void* addr) {
    _address = addr;
    _value   = *((T*)_address);
  }

  /**
   * getValue -> Returns the value of the variable
   */
  double getValue() {
    //if we are dealing with a reference, it might need to refresh itself
    if (_address != NULL)
      _value   = *((T*)_address);
    return (double) _value;
  }
};

}
#endif
