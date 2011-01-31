#include <math.h>
#include <cmath>
#include <ctime>
#include <math.h>       //Needed for erf and erfc on solaris.

#include <stk_expreval/Function.hpp>
#include <stk_expreval/Constants.hpp>


namespace stk {
namespace expreval {

extern "C" {
  typedef double (*CExtern0)();
  typedef double (*CExtern1)(double);
  typedef double (*CExtern2)(double, double);
}


template <>
class CFunction<CExtern0> : public CFunctionBase
{
public:
  typedef CExtern0 Signature;

  explicit CFunction<CExtern0>(Signature function)
    : CFunctionBase(0),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)();
  }

private:
  Signature	m_function;
};


template <>
class CFunction<CExtern1> : public CFunctionBase
{
public:
  typedef CExtern1 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(1),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)(argv[0]);
  }

private:
  Signature	m_function;
};


template <>
class CFunction<CExtern2> : public CFunctionBase
{
public:
  typedef CExtern2 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(2),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)(argv[0], argv[1]);
  }

private:
  Signature	m_function;
};


typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;


extern "C" {
  static double ipart(double x)  {
    double y;
    std::modf(x, &y);
    return y;
  }

  static double fpart(double x)  {
    double y;
    return std::modf(x, &y);
  }

  static double real_rand()  {
    return (double) std::rand() / ((double)(RAND_MAX) + 1.0);
  }

  static double real_srand(double x)  {
    std::srand((int) x);
    return 0.0;
  }

  static double randomize()  {
    std::srand(::time(NULL));
    return 0.0;
  }

  static double deg(double a)  {
    return (180.0 / s_pi) * a;
  }

  static double rad(double a)  {
    return  (s_pi / 180.0) * a;
  }

  static double min(double x, double y) {
    return std::min(x, y);
  }

  static double max(double x, double y) {
    return std::max(x, y);
  }

  static double recttopolr(double x, double y) {
    return std::sqrt((x * x) + (y * y));
  }

  static double recttopola(double x, double y) {
    double tmp = std::atan2(y, x);

    /* Convert to 0.0 to 2 * PI */
    if (tmp < 0.0)
      return tmp + (2.0 * s_pi);
    else
      return tmp;
  }

  static double poltorectx(double r, double theta) {
    return r * std::cos(theta);
  }

  static double poltorecty(double r, double theta) {
    return r * std::sin(theta);
  }
}


CFunctionMap::CFunctionMap() 
{
  (*this)["rand"] = new CFunction0(real_rand);
  (*this)["srand"] = new CFunction1(real_srand);
  (*this)["randomize"] = new CFunction0(randomize);

  (*this)["exp"] = new CFunction1(std::exp);
  (*this)["ln"] = new CFunction1(std::log);
  (*this)["log"] = new CFunction1(std::log);
  (*this)["log10"] = new CFunction1(std::log10);
  (*this)["pow"] = new CFunction2(std::pow);
  (*this)["sqrt"] = new CFunction1(std::sqrt);
  (*this)["erfc"] = new CFunction1(erfc);
  (*this)["erf"] = new CFunction1(erf);

  (*this)["acos"] = new CFunction1(std::acos);
  (*this)["asin"] = new CFunction1(std::asin);
  (*this)["atan"] = new CFunction1(std::atan);
  (*this)["atan2"] = new CFunction2(std::atan2);
  (*this)["ceil"] = new CFunction1(std::ceil);
  (*this)["cos"] = new CFunction1(std::cos);
  (*this)["cosh"] = new CFunction1(std::cosh);
  (*this)["floor"] = new CFunction1(std::floor);
  (*this)["sin"] = new CFunction1(std::sin);
  (*this)["sinh"] = new CFunction1(std::sinh);
  (*this)["tan"] = new CFunction1(std::tan);
  (*this)["tanh"] = new CFunction1(std::tanh);

  (*this)["abs"] = new CFunction1(std::fabs);
  (*this)["fabs"] = new CFunction1(std::fabs);
  (*this)["deg"] = new CFunction1(deg);
  (*this)["mod"] = new CFunction2(std::fmod);
  (*this)["fmod"] = new CFunction2(std::fmod);
  (*this)["ipart"] = new CFunction1(ipart);
  (*this)["fpart"] = new CFunction1(fpart);
  (*this)["max"] = new CFunction2(max);
  (*this)["min"] = new CFunction2(min);
  (*this)["poltorectx"] = new CFunction2(poltorectx);
  (*this)["poltorecty"] = new CFunction2(poltorecty);
  (*this)["rad"] = new CFunction1(rad);
  (*this)["recttopola"] = new CFunction2(recttopola);
  (*this)["recttopolr"] = new CFunction2(recttopolr);
}


CFunctionMap::~CFunctionMap()
{
  for (CFunctionMap::iterator it = begin(); it != end(); ++it)
    delete (*it).second;
}


CFunctionMap &
getCFunctionMap()
{
  static CFunctionMap s_functionMap;

  return s_functionMap;
}

} // namespace expreval
} // namespace stk
