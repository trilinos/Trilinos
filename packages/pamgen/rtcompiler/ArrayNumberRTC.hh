#ifndef _ARRAYNUMBERRTC_H
#define _ARRAYNUMBERRTC_H

#include <iostream>
#include <assert.h>
#include "ValueRTC.hh"
#include "commonRTC.hh"

namespace PG_RuntimeCompiler {

/**
 * A Number object represnts the constants in the code the user gives us. 
 */

template <class T> 
class ArrayNumber: public Value
{
 public:

  /**
   * Constructor -> Trivial
   *
   * @param type  - The type of the number
   * @param value - The value of the number
   */
  ArrayNumber(Type type, const T* value, int size): Value(type, ArrayNumberOT) 
  { 
    assert(size > 0);
    assert(value != NULL);

    _size  = size;
    _value = new T[size];
    for (int i = 0; i < size; ++i) 
      _value[i] = value[i];
  }
  
  ~ArrayNumber() { delete[] _value; }

  /**
   * getValue -> Returns the value of the number
   */
  double getArrayValue(int offset) const {
    assert(offset >= 0 && offset < _size);

    return (double) _value[offset];
  }

  /**
   * setValue -> Changes the value of the number
   *
   * @param value - The new value for the number
   */ 
  void setArrayValue(double value, int offset) {
    assert(offset >= 0 && offset < _size);

    _value[offset] = (T) value;
  }

  int getSize() const {return _size;}

  void print() { 
    std::cout << "{";
    for (int i = 0; i < _size; ++i)
      std::cout << _value[i] << ", ";
    std::cout << "}" << std::endl;
  }

 protected:
  T*  _value; //!< The value of this number
  int _size;
};

}

#endif
