// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_OperatorRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_VariableRTC.hh"
#include "RTC_ObjectRTC.hh"
#include "RTC_ScalarNumberRTC.hh"

#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace std;
using namespace PG_RuntimeCompiler;

bool Operator::ISINIT = false;
map<string, Operator*> Operator::OPERATORS;

/*****************************************************************************/
Operator::Operator(const string& name, bool is_Unary,
                   int precedence, bool isLeftAssoc) : Object(OperatorOT)
/*****************************************************************************/
{
  _name            = name;
  _isUnary         = is_Unary;
  _precedence      = precedence;
  _leftAssociative = isLeftAssoc;
}

/*****************************************************************************/
void Operator::init()
/*****************************************************************************/
{
  if (!ISINIT) {
    static AddO          o1;  OPERATORS[o1.getName()]  = &o1;
    static SubtractO     o2;  OPERATORS[o2.getName()]  = &o2;
    static MultiplyO     o3;  OPERATORS[o3.getName()]  = &o3;
    static DivideO       o4;  OPERATORS[o4.getName()]  = &o4;
    static EqualityO     o5;  OPERATORS[o5.getName()]  = &o5;
    static GreaterO      o6;  OPERATORS[o6.getName()]  = &o6;
    static LessO         o7;  OPERATORS[o7.getName()]  = &o7;
    static GreaterEqualO o8;  OPERATORS[o8.getName()]  = &o8;
    static LessEqualO    o9;  OPERATORS[o9.getName()]  = &o9;
    static AssignmentO   o10; OPERATORS[o10.getName()] = &o10;
    static LogicOrO      o11; OPERATORS[o11.getName()] = &o11;
    static LogicAndO     o12; OPERATORS[o12.getName()] = &o12;
    static InEqualityO   o13; OPERATORS[o13.getName()] = &o13;
    static ModuloO       o14; OPERATORS[o14.getName()] = &o14;
    static ExponentO     o15; OPERATORS[o15.getName()] = &o15;
    static OpenParenO    o16; OPERATORS[o16.getName()] = &o16;
    static CloseParenO   o17; OPERATORS[o17.getName()] = &o17;
    static ArrayInitO    o18; OPERATORS[o18.getName()] = &o18;
    static NegationO     o19; OPERATORS[o19.getName()] = &o19;
    static LogicNotO     o20; OPERATORS[o20.getName()] = &o20;

    ISINIT = true;
  }
}

/*****************************************************************************/
Operator* Operator::getOp(const string& name)
/*****************************************************************************/
{
  if (!ISINIT) init();

  return (OPERATORS.find(name) != OPERATORS.end()) ? OPERATORS[name] : NULL;
}

/*****************************************************************************/
ostream& Operator::operator<<(ostream& os) const
/*****************************************************************************/
{
  os << getName();
  return os;
}

/*****************************************************************************/
void AddO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() + arg2->getValue());
}

/*****************************************************************************/
void SubtractO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() - arg2->getValue());
}

/*****************************************************************************/
void MultiplyO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() * arg2->getValue());
}

/*****************************************************************************/
void DivideO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() / arg2->getValue());
}

/*****************************************************************************/
void EqualityO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() == arg2->getValue());
}

/*****************************************************************************/
void GreaterO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() > arg2->getValue());
}

/*****************************************************************************/
void LessO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() < arg2->getValue());
}

/*****************************************************************************/
void GreaterEqualO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() >= arg2->getValue());
}

/*****************************************************************************/
void LessEqualO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() <= arg2->getValue());
}

/*****************************************************************************/
void AssignmentO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  assert(isAssignable(arg1->getObjectType()));
  arg1->setValue(arg2->getValue());
  store.setValue(arg2->getValue());
}

/*****************************************************************************/
void LogicOrO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() || arg2->getValue());
}

/*****************************************************************************/
void LogicAndO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() && arg2->getValue());
}

/*****************************************************************************/
void InEqualityO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(arg1->getValue() != arg2->getValue());
}

/*****************************************************************************/
void ModuloO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue((int)arg1->getValue() % (int)arg2->getValue());
}

/*****************************************************************************/
void ExponentO::performOp(Value* arg1, Value* arg2, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(pow(arg1->getValue(), arg2->getValue()));
}

/*****************************************************************************/
void NegationO::performOp(Value* arg, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(-arg->getValue());
}

/*****************************************************************************/
void LogicNotO::performOp(Value* arg, ScalarNumber<double>& store)
/*****************************************************************************/
{
  store.setValue(!arg->getValue());
}
