// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _REGISTRAR_H
#define _REGISTRAR_H

#include "RTC_ValueRTC.hh"
#include "RTC_ObjectRTC.hh"
#include "RTC_LineRTC.hh"

#include <map>
#include <list>
#include <string>
#include <iostream>
#include <cassert>

namespace PG_RuntimeCompiler {

/**
 * RTBoundFunc is the class all callable function classes must inherit from.
 */
class RTBoundFunc
{
 public:
  RTBoundFunc(const std::string& name_, int args, bool is_Optimizable = true,
              bool variable_NumArgs = false) {
    _name             = name_;
    _numArgs          = args;
    _isOptimizable    = is_Optimizable;
    _variableNumArgs  = variable_NumArgs;
  }
  virtual ~RTBoundFunc() {}
  virtual double execute(Value**) = 0;

  bool isOptimizable() const {return _isOptimizable;}

  bool variableNumArgs() const {return _variableNumArgs;}

  std::string name() const {return _name;}

  int numArgs() const {return _numArgs;}

  void numArgs(int args) {
    assert(_variableNumArgs);
    _numArgs = args;
  }

 protected:
  std::string _name;
  int         _numArgs; //! In the case of variableNumArgs, this is a minimum
  bool        _isOptimizable;
  bool        _variableNumArgs;
};

/**
 * This class is the super class for all function calls.
 */
class FunctionCall : public Object
{
 public:
  FunctionCall() : Object(FunctionOT) {}

  virtual ~FunctionCall() {}

  virtual double execute() = 0;

  virtual void fillArg(Line* l) = 0;

  virtual bool canGoEarly() const = 0;

  bool hasVariableArgs() const {return _func->variableNumArgs();}
  unsigned getNumArgs() const { return _func->numArgs(); }

protected:
  RTBoundFunc* _func;
};

/**
 * The FixedArgFunctionCall class represents individual function calls in the
 * user's code. This class is for functions with a fixed number of arguments.
 */
class FixedArgFunctionCall : public FunctionCall
{
 public:
  FixedArgFunctionCall(RTBoundFunc* rt);

  ~FixedArgFunctionCall();

  std::ostream& operator<<(std::ostream& os) const;

  double execute();

  void fillArg(Line* l);

  bool canGoEarly() const;

 protected:
  int          _argIndex;
  Line**       _argExpressions;
  Value**      _argValues;
};

/**
 * The VariableArgFunctionCall class represents individual function calls
 * in the user's code. This class is for functions that can have a variable
 * number of arguments, like printf.
 */
class VariableArgFunctionCall : public FunctionCall
{
 public:
  VariableArgFunctionCall(RTBoundFunc* rt);

  ~VariableArgFunctionCall();

  std::ostream& operator<<(std::ostream& os) const;

  double execute();

  void fillArg(Line* l);

  bool canGoEarly() const;

 private:
  std::list<Line*> _argExpressionList;
  Value**          _argValues;
};


/**
 * Registrar is used to contain all the registered functions that can be called
 * by users.
 */
class Registrar
{
 private:
  static std::map<std::string,RTBoundFunc*> FUNCTIONS;
  static bool ISINIT;

 public:
  static void registerFunction(RTBoundFunc*);
  static RTBoundFunc* findByName(const std::string& );
  static void setupStandardFunctions();
  static FunctionCall* generateCall(const std::string &);
};

class Fabs : public RTBoundFunc
{
 public:
  Fabs() : RTBoundFunc("fabs", 1) {}

  double execute(Value**);
};

class Sin : public RTBoundFunc
{
 public:
  Sin() : RTBoundFunc("sin", 1) {}

  double execute(Value**);
};

class Tan : public RTBoundFunc
{
 public:
  Tan() : RTBoundFunc("tan", 1) {}

  double execute(Value**);
};

class Cos : public RTBoundFunc
{
 public:
  Cos() : RTBoundFunc("cos", 1) {}

  double execute(Value**);
};

class Sqrt : public RTBoundFunc
{
 public:
  Sqrt() : RTBoundFunc("sqrt", 1) {}

  double execute(Value**);
};

class Tanh : public RTBoundFunc
{
 public:
  Tanh() : RTBoundFunc("tanh", 1) {}

  double execute(Value**);
};

class Sinh : public RTBoundFunc
{
 public:
  Sinh() : RTBoundFunc("sinh", 1) {}

  double execute(Value**);
};

class Cosh : public RTBoundFunc
{
 public:
  Cosh() : RTBoundFunc("cosh", 1) {}

  double execute(Value**);
};

class Log : public RTBoundFunc
{
 public:
  Log() : RTBoundFunc("log", 1) {}

  double execute(Value**);
};

class Log10 : public RTBoundFunc
{
 public:
  Log10() : RTBoundFunc("log10", 1) {}

  double execute(Value**);
};

class Exp : public RTBoundFunc
{
 public:
  Exp() : RTBoundFunc("exp", 1) {}

  double execute(Value**);
};

class Atan : public RTBoundFunc
{
 public:
  Atan() : RTBoundFunc("atan", 1) {}

  double execute(Value**);
};

class Atantwo : public RTBoundFunc
{
 public:
  Atantwo() : RTBoundFunc("atan2", 2) {}

  double execute(Value**);
};

class Asin : public RTBoundFunc
{
 public:
  Asin() : RTBoundFunc("asin", 1) {}

  double execute(Value**);
};

class Acos : public RTBoundFunc
{
 public:
  Acos() : RTBoundFunc("acos", 1) {}

  double execute(Value**);
};

class Pow : public RTBoundFunc
{
 public:
  Pow() : RTBoundFunc("pow", 2) {}

  double execute(Value**);
};

class Rand : public RTBoundFunc
{
 public:
  Rand() : RTBoundFunc("rand", 0, false) {}

  double execute(Value**);
};

class DRand : public RTBoundFunc
{
 public:
  DRand() : RTBoundFunc("drand", 0, false) {}

  double execute(Value**);
};

class Gamma : public RTBoundFunc
{
 public:
  Gamma() : RTBoundFunc("gamma", 1) {}

  double execute(Value**);
};

class Print : public RTBoundFunc
{
 public:
  Print() : RTBoundFunc("print", 1, false) {}

  double execute(Value**);
};


class Printf : public RTBoundFunc
{
 public:
  Printf() : RTBoundFunc("printf", 1, false, true) {}

  double execute(Value**);
};

class Tester : public RTBoundFunc
{
 public:
  Tester() : RTBoundFunc("arrayFill", 2) {}

  double execute(Value**);
};

class Bessel_J0 : public RTBoundFunc
{
 public:
  Bessel_J0() : RTBoundFunc("j0", 1) {}

  double execute(Value**);
};

class Bessel_J1 : public RTBoundFunc
{
 public:
  Bessel_J1() : RTBoundFunc("j1", 1) {}

  double execute(Value**);
};

class Bessel_I0 : public RTBoundFunc
{
 public:
  Bessel_I0() : RTBoundFunc("i0", 1) {}

  double execute(Value**);
};

class Bessel_I1 : public RTBoundFunc
{
 public:
  Bessel_I1() : RTBoundFunc("i1", 1) {}

  double execute(Value**);
};

class Erf : public RTBoundFunc
{
 public:
  Erf() : RTBoundFunc("erf", 1) {}

  double execute(Value**);
};

class Erfc : public RTBoundFunc
{
 public:
  Erfc() : RTBoundFunc("erfc", 1) {}

  double execute(Value**);
};

class GeneralizedCompleteEllipticIntegral : public RTBoundFunc
{
 public:
  GeneralizedCompleteEllipticIntegral () : RTBoundFunc("gcei", 4) {}

  double execute(Value**);
};

class Readline : public RTBoundFunc
{
 public:
  // read(filename, buffer), returns number of chars read, returns 0 if end of file
  Readline() : RTBoundFunc("readline", 2) {}

  double execute(Value**);
};

class Scanf : public RTBoundFunc
{
 public:
  Scanf() : RTBoundFunc("scanf", 1, false, true) {}

  double execute(Value**);
};

}
#endif
