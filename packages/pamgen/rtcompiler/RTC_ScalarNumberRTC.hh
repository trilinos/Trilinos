// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _SCALARNUMBERRTC_H
#define _SCALARNUMBERRTC_H

#include "RTC_ValueRTC.hh"
#include "RTC_commonRTC.hh"

#include <iostream>
#include <cassert>

namespace PG_RuntimeCompiler {

/**
 * A Number object represnts the constants in the code the user gives us.
 */

template <class T>
class ScalarNumber: public Value
{
 public:

  /**
   * Constructor -> Trivial
   *
   * @param value - The value of the number
   */
  ScalarNumber(T value) : Value(TypeToTypeT<T>::value, ScalarNumberOT) {
    _value = value;
  }

  /**
   * Constructor -> Value set to zero
   */
  ScalarNumber() : Value(TypeToTypeT<T>::value, ScalarNumberOT) { _value = 0; }

  /**
   * getValue -> Returns the value of the number
   */
  double getValue()  {return (double) _value;}

  /**
   * setValue -> Changes the value of the number
   *
   * @param value - The new value for the number
   */
  void setValue(double value) {_value = (T) value;}

  std::ostream& operator<<(std::ostream& os) const {
    //os << "ScalarNumber:" << _value;
    os << _value;
    return os;
  }

 protected:
  T _value; //!< The value of this number
};

}

#endif
