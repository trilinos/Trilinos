// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _OPERATORRTC_H
#define _OPERATORRTC_H

#include "RTC_commonRTC.hh"
#include "RTC_ObjectRTC.hh"
#include "RTC_ScalarNumberRTC.hh"

#include <string>
#include <map>
#include <cassert>


namespace PG_RuntimeCompiler {

/**
 * A Operator object represnts the operators in the code the user gives us.
 */

class Value;

class Operator: public Object
{
 public:

  /**
   * Constructor -> Trivial
   *
   * @param op - An enum specifying which operator this is
   */
  Operator(const std::string& name, bool isUnary,
           int precedence, bool isLeftAssoc);

  /**
   * Destructor -> The destructor is a no-op
   */
  virtual ~Operator() {}

  /**
   * getName() -> This method returns the string name of the
   *              Operator(ex: "Add","Subtract",...)
   */
  std::string getName() const { return _name;}

  int getPrecedence() const {return _precedence;}

  bool isLeftAssociative() const {return _leftAssociative;}

  /**
   * doBinaryOp -> This method performs an operation on arg1 and arg2 and
   *               stores the result in store.
   *
   * @param arg1  - The first  argument to the operation
   * @param arg2  - The second argument to the operation
   * @param store - The memory where the result is stored
   */
  virtual void performOp(Value* arg1, Value* arg2,
                         ScalarNumber<double>& store) = 0;

  virtual void performOp(Value* arg, ScalarNumber<double>& store) = 0;

  std::ostream& operator<<(std::ostream& os) const;

  bool isUnary() const {return _isUnary;}

  static Operator* getOp(const std::string& name);

  static void init();
 protected:

  bool        _isUnary;
  std::string _name; //!< A string representation of the operator
  int         _precedence;
  bool        _leftAssociative;

 private:
  static std::map<std::string, Operator*> OPERATORS;
  static bool ISINIT;
};

class AddO : public Operator
{
 public:
  AddO() : Operator("+", false, 6, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class SubtractO : public Operator
{
 public:
  SubtractO() : Operator("-", false, 6, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class MultiplyO : public Operator
{
 public:
  MultiplyO() : Operator("*", false, 7, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class DivideO : public Operator
{
 public:
  DivideO() : Operator("/", false, 7, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class EqualityO : public Operator
{
 public:
  EqualityO() : Operator("==", false, 3, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class GreaterO : public Operator
{
 public:
  GreaterO() : Operator(">", false, 4, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class LessO : public Operator
{
 public:
  LessO() : Operator("<", false, 4, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class GreaterEqualO : public Operator
{
 public:
  GreaterEqualO() : Operator(">=", false, 4, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class LessEqualO : public Operator
{
 public:
  LessEqualO() : Operator("<=", false, 4, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class AssignmentO : public Operator
{
 public:
  AssignmentO() : Operator("=", false, 0, false) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class LogicOrO : public Operator
{
 public:
  LogicOrO() : Operator("||", false, 1, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class LogicAndO : public Operator
{
 public:
  LogicAndO() : Operator("&&", false, 2, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class InEqualityO : public Operator
{
 public:
  InEqualityO() : Operator("!=", false, 3, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class ModuloO : public Operator
{
 public:
  ModuloO() : Operator("%", false, 7, true) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class ExponentO : public Operator
{
 public:
  ExponentO() : Operator("^", false, 8, false) {}

  void performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store);
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

//Important - OpenParen, CloseParen, and ArrayInit are structural operators,
//not arithmetic operators, therefore they cannot be performed.
class OpenParenO : public Operator
{
 public:
  OpenParenO() : Operator("(", false, 0, true) {}

  void performOp(Value* /* arg1 */, Value* /* arg2 */, ScalarNumber<double>& /* store */) {assert(0);}
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class CloseParenO : public Operator
{
 public:
  CloseParenO() : Operator(")", false, 0, true) {}

  void performOp(Value* /* arg1 */, Value* /* arg2 */, ScalarNumber<double>& /* store */){assert(0);}
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class ArrayInitO : public Operator
{
 public:
  ArrayInitO() : Operator("init", false, 0, true) {}

  void performOp(Value* /* arg1 */, Value* /* arg2 */, ScalarNumber<double>& /* store */){assert(0);}
  void performOp(Value* /* arg */, ScalarNumber<double>& /* store */) {assert(0);}
};

class NegationO : public Operator
{
 public:
  NegationO() : Operator("_", true, 11, false) {}

  void performOp(Value* /* arg1 */, Value* /* arg2 */, ScalarNumber<double>& /* store */) {assert(0);}
  void performOp(Value* arg, ScalarNumber<double>& store);
};

class LogicNotO : public Operator
{
 public:
  LogicNotO() : Operator("!", true, 11, false) {}

  void performOp(Value* /* arg1 */, Value* /* arg2 */, ScalarNumber<double>& /* store */) {assert(0);}
  void performOp(Value* arg, ScalarNumber<double>& store);
};

}

#endif
