#include "Pamgen_config.h"
#include "RTC_LineRTC.hh"
#include "RTC_RegistrarRTC.hh"
#include "RTC_ArrayVarRTC.hh"
#include "RTC_Bessel_I_RTC.hh"
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>

#ifdef HAVE_PAMGEN_BOOST
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#endif

using namespace std;
using namespace PG_RuntimeCompiler;

map<string,RTBoundFunc*> Registrar::FUNCTIONS;
bool Registrar::ISINIT = false;

/*****************************************************************************/
FixedArgFunctionCall::FixedArgFunctionCall(RTBoundFunc* rt) : FunctionCall()
/*****************************************************************************/
{
  assert(!rt->variableNumArgs());
  _func     = rt;
  _argIndex = 0;

  if (_func->numArgs() != 0) {
    _argExpressions = new Line*[_func->numArgs()];
    _argValues      = new Value*[_func->numArgs()];
    for (int i = 0; i < _func->numArgs(); ++i) {
      _argExpressions[i] = NULL;
      _argValues[i]      = NULL;
    }
  }
  else {
    _argExpressions = NULL;
    _argValues      = NULL;
  }
}

/*****************************************************************************/
FixedArgFunctionCall::~FixedArgFunctionCall()
/*****************************************************************************/
{
  if (_argExpressions != NULL) {
    assert(_argValues != NULL);

    for (int i = 0; i < _func->numArgs(); ++i)
      delete _argExpressions[i];

    delete[] _argExpressions;
    delete[] _argValues;
  }
}

/*****************************************************************************/
ostream& FixedArgFunctionCall::operator<<(ostream& os) const
/*****************************************************************************/
{
  os << _func->name() << "(";
  for (int i = 0; i < _func->numArgs(); ++i) {
    os << *(_argExpressions[i]);
    if (i < _func->numArgs() - 1) os << ", ";
  }
  os << ")";

  return os;
}

/*****************************************************************************/
double FixedArgFunctionCall::execute()
/*****************************************************************************/
{
  //Get the values of the argument expressions
  for (int i = 0; i < _func->numArgs(); ++i) {
    assert(_argExpressions[i] != NULL);
    _argValues[i] = _argExpressions[i]->execute();
  }

  return _func->execute(_argValues);
}

/*****************************************************************************/
void FixedArgFunctionCall::fillArg(Line* line)
/*****************************************************************************/
{
  assert(line != NULL);
  assert(_argIndex < _func->numArgs());
  _argExpressions[_argIndex++] = line;
}

/*****************************************************************************/
bool FixedArgFunctionCall::canGoEarly() const
/*****************************************************************************/
{
  if (!_func->isOptimizable())
    return false;

  for (int i = 0; i < _func->numArgs(); ++i) {
    assert(_argExpressions[i] != NULL);
    if (!_argExpressions[i]->can_go())
      return false;
  }
  return true;
}

/*****************************************************************************/
VariableArgFunctionCall::VariableArgFunctionCall(RTBoundFunc* rt)
  : FunctionCall()
/*****************************************************************************/
{
  assert(rt->variableNumArgs());

  _func      = rt;
  _argValues = NULL;
}

/*****************************************************************************/
VariableArgFunctionCall::~VariableArgFunctionCall()
/*****************************************************************************/
{
  if (_argExpressionList.size() > 0) {

    list<Line*>::iterator itr = _argExpressionList.begin();
    for ( ; itr != _argExpressionList.end(); ++itr)
      delete (*itr);
    if(_argValues != 0)
      delete[] _argValues;
  }
}

/*****************************************************************************/
ostream& VariableArgFunctionCall::operator<<(ostream& os) const
/*****************************************************************************/
{
  os << _func->name() << "(";

  list<Line*>::const_iterator itr = _argExpressionList.begin();
  for ( ; itr != _argExpressionList.end(); ++itr)
    os << *(*itr) << ", ";

  os << ")";

  return os;
}

/*****************************************************************************/
double VariableArgFunctionCall::execute()
/*****************************************************************************/
{
  //Get the values of the argument expressions
  if (_argValues == NULL && _argExpressionList.size() > 0) {
    _argValues = new Value*[_argExpressionList.size()];
    _func->numArgs(_argExpressionList.size());
  }

  list<Line*>::iterator itr = _argExpressionList.begin();
  for (int i = 0; itr != _argExpressionList.end(); ++itr, ++i)
    _argValues[i] = (*itr)->execute();

  return _func->execute(_argValues);
}

/*****************************************************************************/
void VariableArgFunctionCall::fillArg(Line* line)
/*****************************************************************************/
{
  _argExpressionList.push_back(line);
}

/*****************************************************************************/
bool VariableArgFunctionCall::canGoEarly() const
/*****************************************************************************/
{
  if (!_func->isOptimizable())
    return false;

  list<Line*>::const_iterator itr = _argExpressionList.begin();
  for ( ; itr != _argExpressionList.end(); ++itr) {
    if (!(*itr)->can_go())
      return false;
  }
  return true;
}


/*****************************************************************************/
void Registrar::registerFunction(RTBoundFunc* func)
/*****************************************************************************/
{
  FUNCTIONS[func->name()] = func;
}

/*****************************************************************************/
RTBoundFunc* Registrar::findByName(const string& str)
/*****************************************************************************/
{
  if (!ISINIT) setupStandardFunctions();

  return (FUNCTIONS.find(str) != FUNCTIONS.end()) ? FUNCTIONS[str] : NULL;
}

/*****************************************************************************/
FunctionCall* Registrar::generateCall(const string& str)
/*****************************************************************************/
{
  RTBoundFunc* rt = findByName(str);
  if (rt == NULL)
    return NULL;

  if (rt->variableNumArgs())
    return new VariableArgFunctionCall(rt);
  else
    return new FixedArgFunctionCall(rt);
}

/*****************************************************************************/
void Registrar::setupStandardFunctions()
/*****************************************************************************/
{
  if (!ISINIT) {
    static Fabs      f1;  FUNCTIONS[f1.name()]  = &f1;
    static Sin       f2;  FUNCTIONS[f2.name()]  = &f2;
    static Cos       f3;  FUNCTIONS[f3.name()]  = &f3;
    static Tan       f4;  FUNCTIONS[f4.name()]  = &f4;
    static Sqrt      f5;  FUNCTIONS[f5.name()]  = &f5;
    static Tanh      f6;  FUNCTIONS[f6.name()]  = &f6;
    static Sinh      f7;  FUNCTIONS[f7.name()]  = &f7;
    static Cosh      f8;  FUNCTIONS[f8.name()]  = &f8;
    static Log       f9;  FUNCTIONS[f9.name()]  = &f9;
    static Log10     f10; FUNCTIONS[f10.name()] = &f10;
    static Exp       f11; FUNCTIONS[f11.name()] = &f11;
    static Asin      f12; FUNCTIONS[f12.name()] = &f12;
    static Atan      f13; FUNCTIONS[f13.name()] = &f13;
    static Atantwo   f14; FUNCTIONS[f14.name()] = &f14;
    static Acos      f15; FUNCTIONS[f15.name()] = &f15;
    static Rand      f16; FUNCTIONS[f16.name()] = &f16;
    static DRand     f17; FUNCTIONS[f17.name()] = &f17;
    static Pow       f18; FUNCTIONS[f18.name()] = &f18;
    static Printf    f19; FUNCTIONS[f19.name()] = &f19;
    static Tester    f20; FUNCTIONS[f20.name()] = &f20;
    static Bessel_J0 f21; FUNCTIONS[f21.name()] = &f21;
    static Bessel_J1 f22; FUNCTIONS[f22.name()] = &f22;
    static Bessel_I0 f23; FUNCTIONS[f23.name()] = &f23;
    static Bessel_I1 f24; FUNCTIONS[f24.name()] = &f24;
    static Erf       f25; FUNCTIONS[f25.name()] = &f25;
    static Erfc      f26; FUNCTIONS[f26.name()] = &f26;
    static Gamma     f27; FUNCTIONS[f27.name()] = &f27;
    static Print     f28; FUNCTIONS[f28.name()] = &f28;
    static GeneralizedCompleteEllipticIntegral      f29; FUNCTIONS[f29.name()] = &f29;
    static Readline  f30; FUNCTIONS[f30.name()] = &f30;
    static Scanf     f31; FUNCTIONS[f31.name()] = &f31;
    ISINIT = true;
  }
}

/*****************************************************************************/
double Fabs::execute(Value** args)
/*****************************************************************************/
{
  return fabs(args[0]->getValue());
}

/*****************************************************************************/
double Sin::execute(Value** args)
/*****************************************************************************/
{
  return sin(args[0]->getValue());
}

/*****************************************************************************/
double Cos::execute(Value** args)
/*****************************************************************************/
{
  return cos(args[0]->getValue());
}

/*****************************************************************************/
double Tan::execute(Value** args)
/*****************************************************************************/
{
  return tan(args[0]->getValue());
}

/*****************************************************************************/
double Sqrt::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getValue() >= 0); //cant take sqrt of a negative
  return sqrt(args[0]->getValue());
}

/*****************************************************************************/
double Tanh::execute(Value** args)
/*****************************************************************************/
{
  return tanh(args[0]->getValue());
}

/*****************************************************************************/
double Sinh::execute(Value** args)
/*****************************************************************************/
{
  return sinh(args[0]->getValue());
}

/*****************************************************************************/
double Cosh::execute(Value** args)
/*****************************************************************************/
{
  return cosh(args[0]->getValue());
}

/*****************************************************************************/
double Log::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getValue() >= 0); //cant take log of a negative
  return log(args[0]->getValue());
}

/*****************************************************************************/
double Log10::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getValue() >= 0); //cant take log of a negative
  return log10(args[0]->getValue());
}

/*****************************************************************************/
double Exp::execute(Value** args)
/*****************************************************************************/
{
  return exp(args[0]->getValue());
}

/*****************************************************************************/
double Atan::execute(Value** args)
/*****************************************************************************/
{
  return atan(args[0]->getValue());
}

/*****************************************************************************/
double Atantwo::execute(Value** args)
/*****************************************************************************/
{
  return atan2(args[0]->getValue(), args[1]->getValue());
}

/*****************************************************************************/
double Asin::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getValue() >= -1 && args[0]->getValue() <= 1);
  return asin(args[0]->getValue());
}

/*****************************************************************************/
double Acos::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getValue() >= -1 && args[0]->getValue() <= 1);
  return acos(args[0]->getValue());
}

/*****************************************************************************/
double Pow::execute(Value** args)
/*****************************************************************************/
{
  return pow(args[0]->getValue(), args[1]->getValue());
}

/*****************************************************************************/
double Rand::execute(Value** args)
/*****************************************************************************/
{
  assert(args == NULL);
  return rand();
}

/*****************************************************************************/
double DRand::execute(Value** args)
/*****************************************************************************/
{
  assert(args == NULL);
  //returns random between 0 and 1
  return ((double)rand()) / ((double) RAND_MAX);;
}

/*****************************************************************************/
double Gamma::execute(Value** args)
/*****************************************************************************/
{
#if _MSC_VER
  return 0;
#else
  return tgamma(args[0]->getValue());
#endif
}

/*****************************************************************************/
double Print::execute(Value** args)
/*****************************************************************************/
{
  cout << *(args[0]) << endl;
  return 0;
}

/*****************************************************************************/
double Printf::execute(Value** args)
/*****************************************************************************/
{
  //first argument to printf should be a character array
  assert(args[0]->getObjectType() == ArrayNumberOT);
  assert(args[0]->getType() == CharT);

  int numVarArgsUsed = 0;

  for (long i = 0; i < args[0]->getSize(); ++i) {
    if ( ((char) args[0]->getArrayValue(i)) == '%') {
      assert(numVarArgsUsed+1 < _numArgs);
      cout << *(args[numVarArgsUsed+1]);
      numVarArgsUsed++;
    }
    else if ( ((char) args[0]->getArrayValue(i)) == '\"') {
      // Do not print
    }
    else if ( ((char) args[0]->getArrayValue(i)) == '\\' &&
              ((char) args[0]->getArrayValue(i+1)) == '\"') {
      // Do not print
    }
    else {
      cout << ((char) args[0]->getArrayValue(i));
    }
  }
  cout << endl;

  return args[0]->getSize();
}

/*****************************************************************************/
double Tester::execute(Value** args)
/*****************************************************************************/
{
  assert(args[0]->getObjectType() == ArrayVarOT);

  //fills the array in args[0] with the value of args[1]
  for (long i = 0; i < args[0]->getSize(); ++i)
    args[0]->setArrayValue(args[1]->getValue(), i);

  return args[1]->getValue();
}

/*****************************************************************************/
double Bessel_J0::execute(Value** args)
/*****************************************************************************/
{
  return j0(args[0]->getValue());
}

/*****************************************************************************/
double Bessel_J1::execute(Value** args)
/*****************************************************************************/
{
  return j1(args[0]->getValue());
}


/*****************************************************************************/
double Bessel_I0::execute(Value** args)
/*****************************************************************************/
{
  return PAMGEN_NEVADA::Bessel_I0(args[0]->getValue());
}

/*****************************************************************************/
double Bessel_I1::execute(Value** args)
/*****************************************************************************/
{
  return PAMGEN_NEVADA::Bessel_I1(args[0]->getValue());
}

/*****************************************************************************/
double Erf::execute(Value** args)
/*****************************************************************************/
{
#ifdef _MSC_VER
  return 0;
#else
  return erf(args[0]->getValue());
#endif
}

/*****************************************************************************/
double Erfc::execute(Value** args)
/*****************************************************************************/
{
#ifdef _MSC_VER
  return 0;
#else
  return erfc(args[0]->getValue());
#endif
}

/*****************************************************************************/
double GeneralizedCompleteEllipticIntegral::execute(Value** args)
/*****************************************************************************/
{
#ifdef HAVE_PAMGEN_BOOST
  // Generalized complete elliptic integral, implemented using Carlson's functions
  // See Derby and Olbert, "Cylindrical magnets and ideal solenoids," Am. J. Phys. 78, pp. 229-235, 2010.
  // C(kc,p,c,s) = c *R_F(0,kc^2,1) + 1/3 * (s-p*c) * R_J(0,kc^2,1,p)
  double kc2 = args[0]->getValue() * args[0]->getValue();
  return args[2]->getValue() * boost::math::ellint_rf<double,double,double>(0.0,kc2,1.0) +
    (1.0/3.0)*(args[3]->getValue() - args[1]->getValue()*args[2]->getValue())* boost::math::ellint_rj<double,double,double,double>(0.0,kc2,1.0,args[1]->getValue());
#else
  return 0;
#endif
}

/*****************************************************************************/
double Readline::execute(Value** args)
/*****************************************************************************/
{
  //arguments to read should be a character arrays
  assert(args[0]->getObjectType() == ArrayNumberOT);
  assert(args[0]->getType() == CharT);

  assert(args[1]->getObjectType() == ArrayVarOT);
  assert(args[1]->getType() == CharT);

  string filename = "";
  for (long i = 0; i < args[0]->getSize(); ++i) {
    char c = (char) args[0]->getArrayValue(i);
    if (c != '"') {
      filename += c;
    }
  }

  static string curr_filename("");
  static ifstream curr_file;

  if (curr_filename != filename || curr_file.eof()) {
    curr_filename = filename;
    if (curr_file.is_open()) {
      curr_file.close();
    }
    curr_file.open(filename);
  }

  string line;
  if (getline(curr_file, line)) {
     assert(args[1]->getSize() > long(line.size())); // buffer needs to be big enough
    for (size_t i = 0; i < line.size(); ++i) {
      args[1]->setArrayValue(line[i], i);
    }
    args[1]->setArrayValue(0, line.size()); // Null terminate

    return line.size();
  }
  else {
    return -1;
  }
}

/*****************************************************************************/
double Scanf::execute(Value** args)
/*****************************************************************************/
{
  //first argument to scanf should be a character array
  assert(args[0]->getObjectType() == ArrayNumberOT || args[0]->getObjectType() == ArrayVarOT);
  assert(args[0]->getType() == CharT);
  assert(_numArgs < 22);

  string format = "";

  for (int i = 1; i < _numArgs; ++i) {
    format += "%lf";
  }

  string data = "";
  for (long i = 0; i < args[0]->getSize(); ++i) {
    data += (char) args[0]->getArrayValue(i);
  }

  double d[20];

  sscanf(data.c_str(), format.c_str(), &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[9],
         &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19]);

  for (int i = 1; i < _numArgs; ++i) {
    args[i]->setValue(d[i-1]);
  }

  return _numArgs - 1;
}
