// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_FunctionRTC.hh"
#include "RTC_BlockRTC.hh"
#include "RTC_NormalBlockRTC.hh"
#include "RTC_VariableRTC.hh"
#include "RTC_commonRTC.hh"
#include "RTC_ArrayVarRTC.hh"
#include "RTC_ScalarVarRTC.hh"
#include "RTC_TokenizerRTC.hh"

#include <string>
#include <list>
#include <stack>
#include <cmath>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
Function::Function(const unsigned varCount,const string & name):
/*****************************************************************************/
  _name(name),
  _errors(""),
  _mainBlock(NULL),
  _vars(),
  _compiled(false)
{
  _vars.reserve(varCount);
}

/*****************************************************************************/
bool Function::addVar(const string& type, const string& name)
/*****************************************************************************/
{
  //check to make sure the name is valid
  if (name == "") {
    _errors += "Illegal name given for argument\n";
    return false;
  }

  for (unsigned int i = 0; i < name.length(); ++i) {
    if (!isValidVariableChar(name[i])) {
      _errors += "Illegal name given for argument. Make sure you did not include whitespace in the strings passed to addVar.\n";
      return false;
    }
  }

  //create a new variable
  Variable* ptr;
  if (type == "int")
    ptr = new ScalarVar<int>(name);
  else if (type == "long")
    ptr = new ScalarVar<long>(name);
  else if (type == "char")
    ptr = new ScalarVar<char>(name);
  else if (type == "float")
    ptr = new ScalarVar<float>(name);
  else if (type == "double")
    ptr = new ScalarVar<double>(name);
  else if (type == "int[]")
    ptr = new ArrayVar<int>(name);
  else if (type == "long[]")
    ptr = new ArrayVar<long>(name);
  else if (type == "char[]")
    ptr = new ArrayVar<char>(name);
  else if (type == "float[]")
    ptr = new ArrayVar<float>(name);
  else if (type == "double[]")
    ptr = new ArrayVar<double>(name);
  else {
    _errors += "Illegal type provided for argument. Make sure you did not have whitespace in the strings passed to addVar.\n";
    return false;
  }

  //add the variable to the array of argument variables
  _vars.push_back(pair<Variable*,bool>(ptr, false));

  //Mark that it will be initialized by the time the function is executed.
  //We know this because we will crash if any argument has not been filled
  //before the program is executed.
  ptr->init();

  return true;
}

/*****************************************************************************/
bool Function::addBody(const string& body)
/*****************************************************************************/
{
  Tokenizer tokens(body, _errors);

  if (_errors == "") {
    map<string, Variable*> temp; //a map of variables that are in scope

    //add arguments to map of available variables
    for (unsigned int i = 0; i < _vars.size(); ++i)
      temp[_vars[i].first->getName()] = _vars[i].first;

    //Note: additional errors can occur here
    if (_mainBlock != NULL) {
      delete _mainBlock;
    }
    _mainBlock = new NormalBlock(temp, tokens, _errors);
  }

  _compiled = (_errors == "");

  return _compiled;
}

/*****************************************************************************/
void Function::commonVarFill(unsigned index)
/*****************************************************************************/
{
  if (index >= _vars.size()) {
    _errors += "Index is too large.\n";
    return;
  }

  //JGF: commenting this out makes re-fill legal. So, if you fill once, you're
  //free to reuse that fill for N executions, but you are also free to re-fill.
  //if (_vars[index].second) {
  //  _errors += "The variable " + _vars[index].first->getName() +
  //   " was already filled.\n";
  //}

  _vars[index].second = true;
}

/*****************************************************************************/
bool Function::varValueFill(unsigned int index, double value)
/*****************************************************************************/
{
  commonVarFill(index);

  if (_vars[index].first->getObjectType() != ScalarVarOT) {
    _errors += "Must use the varAddrFill method for non-scalar variables.\n";
    return false;
  }

  _vars[index].first->setValue(value);
  return true;
}

/*****************************************************************************/
bool Function::execute()
/*****************************************************************************/
{
  //if addBody did not succeed, we cannot execute
  if (!_compiled) {
    _errors += "Cannot run. Program has not compiled successfully.\n";
    return false;
  }

  //make sure all arguments were filled
  for (unsigned int i = 0; i < _vars.size(); ++i) {
    if (!_vars[i].second) {
      _errors += "The argument " + _vars[i].first->getName() + " was not set.\n";
      return false;
    }
  }

  //execute the function
  _mainBlock->execute();

  //JGF: It is a little unsafe to preserve filled values when the filled value is an
  //address (memory may have been reclaimed).  However, it can also be convenient!
  //mark all variables as unfilled to force refilling before next execution
  //for (unsigned int i = 0; i < _vars.size(); ++i) {
  //  _vars[i].second = false;
  //}

  return true;
}

/*****************************************************************************/
Function::~Function()
/*****************************************************************************/
{
  if (_mainBlock != NULL) {
    delete _mainBlock;
  }

  for (unsigned int i = 0; i < _vars.size(); ++i) {
    delete _vars[i].first;
  }

  _vars.clear();
}

/*****************************************************************************/
void Function::checkType(unsigned int index, int size, int* /* addr */, string& errs)
/*****************************************************************************/
{
  CHECKARGERR((_vars[index].first->getType() != IntT  ||
               (_vars[index].first->getObjectType() == ScalarVarOT && size != 0) ||
               (_vars[index].first->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}

/*****************************************************************************/
void Function::checkType(unsigned int index, int size, long* /* addr */, string& errs)
/*****************************************************************************/
{
  CHECKARGERR((_vars[index].first->getType() != LongT  ||
               (_vars[index].first->getObjectType() == ScalarVarOT && size != 0) ||
               (_vars[index].first->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}

/*****************************************************************************/
void Function::checkType(unsigned int index, int size, float* /* addr */, string& errs)
/*****************************************************************************/
{
  CHECKARGERR((_vars[index].first->getType() != FloatT  ||
               (_vars[index].first->getObjectType() == ScalarVarOT && size != 0) ||
               (_vars[index].first->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}

/*****************************************************************************/
void Function::checkType(unsigned int index, int size, double* /* addr */, string& errs)
/*****************************************************************************/
{
  CHECKARGERR((_vars[index].first->getType() != DoubleT  ||
               (_vars[index].first->getObjectType() == ScalarVarOT && size != 0) ||
               (_vars[index].first->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}

/*****************************************************************************/
void Function::checkType(unsigned int index, int size, char* /* addr */,string& errs)
/*****************************************************************************/
{
  CHECKARGERR((_vars[index].first->getType() != CharT  ||
               (_vars[index].first->getObjectType() == ScalarVarOT && size != 0) ||
               (_vars[index].first->getObjectType() == ArrayVarOT  && size == 0)),
              type_err(index))
}
