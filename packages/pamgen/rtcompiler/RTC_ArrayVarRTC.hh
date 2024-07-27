// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _ARRAYRTC_H
#define _ARRAYRTC_H

#include "RTC_VariableRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_ExecutableRTC.hh"

#include <string>
#include <iostream>

namespace PG_RuntimeCompiler {

/**
 * ArrayVar objects represent variables that are arrays.
 */
template <class T>
class ArrayVar : public Variable
{
 public:

  /**
   * Constructor -> Constructs the super class, initializes instance variables
   *
   * @param name - The name of the variable
   * @param size - The number of elements in the array
   */
  ArrayVar(const std::string& name, long size = 0)
    : Variable(name, TypeToTypeT<T>::value, ArrayVarOT)
  {
    _size    = size;
    _sizeExp = NULL;
    _values  = NULL;
    _isArg   = true; // _values is null
  }

  /**
   * Constructor -> Constructs the super class, initializes instance variables
   *
   * @param name    - The name of the variable
   * @param sizePtr - An expression that, when evaluated, will be
   *                  the array's size
   */
  ArrayVar(const std::string& name, Executable* sizePtr)
    : Variable(name, TypeToTypeT<T>::value, ArrayVarOT)
  {
    _size = 0;
    _sizeExp = sizePtr;
    _values  = NULL;
    _isArg   = true; // _values is null
  }

  /**
   * Destructor -> Delete the size expression if it is not null. We also delete
   *               the values if we are not dealing with user provided
   *               argument.
   */
  virtual ~ArrayVar()
  {
    if (_sizeExp != NULL)
      delete _sizeExp;
    if (!_isArg)
      delete[] _values;
  }

  /**
   * getValue -> Returns the value of the array at a certain index
   *
   * @param offset - The index we want the value of
   */
  double getArrayValue(long offset) const
  {
    assert(_values != NULL);
    if (offset >= _size || offset < 0) {
      std::cout << "Index: " << offset << " is out of bounds on array: "
                << _name << std::endl;
      return 0;
    }
    return (double)(_values[offset]);
  }

  /**
   * setValue -> Sets the value of the array at a certain index
   *
   * @param value  - The value we are going to set to
   * @param offset - The location in the array being changed
   */
  void setArrayValue(double value, long offset)
  {
    assert(_values != NULL);
    if (offset >= _size || offset < 0) {
      std::cout << "Went out of bounds on array: " << _name << std::endl;
      return;
    }
    _values[offset] = (T) value;
  }

  /**
   * getSize -> Returns the size of the array
   */
  long getSize() const {return _size;}

  /**
   * setSize -> Sets the size of the array
   *
   * @param size - The new size of the array
   */
  void setSize(long size) { _size = size; }

  /**
   * evaluateSizeExpr -> Evaluates the size expression to get the array size,
   *                     then allocates a corresponding number of values.
   */
  void evaluateSizeExpr()
  {
    assert(_sizeExp != NULL);
    _size = (int) _sizeExp->execute()->getValue();
    _isArg = false;

    if (_values != NULL)
      delete[] _values;
    _values = new T[_size];
  }

  /**
   * setAddress -> Only called if the array is a function argument. We make the
   *               array point to a user provided address.
   *
   * @param addr - The address of the array's values
   */
  void setAddress(void* addr)
  {
    _isArg  = true;
    _values = (T*) addr;
  }

  std::ostream& operator<<(std::ostream& os) const
  {
    if (_values != NULL) {
      os << "ArrayVar:" << _name << "{";
      for (int i = 0; i < _size; ++i)
        os << _values[i] << ", ";
      os << "}";
    }
    else {
      os << "ArrayVar:" << _name;
    }
    return os;
  }

 private:

  long _size; //!< The length of the array

  T* _values; //!< The array of values

  bool _isArg; /**!< Tells us if the array is a user argument or declared in
                *    in the user defined function. This boolean will affect how
                *    the array memory is cleaned up. We would not want to call
                *    delete on an address the user provided to us as we have no
                *    way of knowing if it points to heap memory.
                */

  Executable* _sizeExp; //!< When evaluated, this will be the size of Array
};

}

#endif
