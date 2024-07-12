// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _ARRAYINDEXRTC_H
#define _ARRAYINDEXRTC_H

#include "RTC_VariableRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_ExecutableRTC.hh"
#include "RTC_ValueRTC.hh"

#include <iostream>
#include <cassert>

namespace PG_RuntimeCompiler {

/**
 * ArrayIndex objects represent the value of an array at a certain index.
 * Anywhere the user get the index of an array,an ArrayIndex object is created.
 * ArrayIndex extends Variable.
 */


class ArrayIndex : public Value {
 public:

  /**
   * Contructor -> Constructs super class, initializes instance variables. Note
   *               this class will take responsibility for deleting index.
   *
   * @param parent - The Array object that is being indexed
   * @param index  - The expression that, when evaluated, will be the index
   */
  ArrayIndex(Variable* parent, Executable* index)
    : Value(parent->getType(), ArrayIndexOT)
  {
    assert(parent->getObjectType() == ArrayVarOT);
    _parent    = parent;
    _indexExpr = index;
  }

  /**
   * Destructor -> Deletes the index expression
   */
  ~ArrayIndex() {delete _indexExpr;}

  /**
   * getValue -> This method evaluates _indexExpr to get the index and then
   *             returns the indexth element of the parent array.
   */
  double getValue() {
    return _parent->getArrayValue((long)_indexExpr->execute()->getValue());
  }

  /**
   * setValue -> This method evaluates _indexExpr to get the index and then
   *             sets the indexth element of the parent array equal to value.
   *
   * @param value - The value we are setting _parent[_indexExpr] equal to
   */
  void setValue(double value) {
    _parent->setArrayValue(value, (long)_indexExpr->execute()->getValue());
  }

  std::ostream& operator<<(std::ostream& os) const {
    //os << "ArrayIndex:" << _parent->getName() << "[" << *_indexExpr << "]";
    os << _parent->getArrayValue((long)_indexExpr->execute()->getValue());
    return os;
  }

 private:

  Variable* _parent;  //!< The Array variable that is being indexed

  Executable* _indexExpr; //!< The expr that, when evaluated, will be the index
};

}
#endif
