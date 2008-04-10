#ifndef _SCALARNUMBERRTC_H
#define _SCALARNUMBERRTC_H

#include <iostream>
#include <assert.h>
#include "ValueRTC.hh"
#include "commonRTC.hh"

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
   * @param type  - The type of the number
   * @param value - The value of the number
   */
  ScalarNumber(Type type, T value) : Value(type, ScalarNumberOT) { 
    _value = value;
  }
  
  /**
   * Constructor -> Value set to zero
   */
  ScalarNumber() : Value(DoubleT, ScalarNumberOT) { _value = 0; }

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

  void print() { std::cout << _value; }
 protected:
  T _value; //!< The value of this number
};

}

#endif
