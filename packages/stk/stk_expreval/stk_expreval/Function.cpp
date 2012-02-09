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
  typedef double (*CExtern3)(double, double, double);
}

static int sRandomRangeHighValue = 3191613;
static int sRandomRangeLowValue  = 1739623;

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

template <>
class CFunction<CExtern3> : public CFunctionBase
{
public:
  typedef CExtern3 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(3),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)(argv[0], argv[1], argv[2]);
  }

private:
  Signature	m_function;
};

typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;
typedef CFunction<CExtern3> CFunction3;


extern "C" {
  /// extract signed integral value from floating-point number
  static double ipart(double x)  {
    double y;
    std::modf(x, &y);
    return y;
  }

  /// Extract fractional value from floating-point number
  static double fpart(double x)  {
    double y;
    return std::modf(x, &y);
  }

  /// Interface to the pseudo-random number generator function rand
  /// provided by ANSI C math library.
  static double real_rand() {
    return (double) std::rand() / ((double)(RAND_MAX) + 1.0);
  }

  /// Sets x as the "seed". Interface to the srand function provided by the
  /// ANSI C math library.
  static double real_srand(double x) {
    std::srand(static_cast<int>(x));
    return 0.0;
  }

  /// Sets the current time as the "seed" to randomize the next call to real_rand.
  static double randomize() {
    std::srand(::time(NULL));
    return 0.0;
  }

  /// Sets x as the "seed" for the pseudo-random number generator.
  static double random_seed(double x) {
    int y = static_cast<int>(x);
    sRandomRangeHighValue =  y;
    sRandomRangeLowValue  = ~y;
    return 0.0;
  }

  /// Non-platform specific (pseudo) random number generator.
  static double random() {
    sRandomRangeHighValue = (sRandomRangeHighValue<<8) + (sRandomRangeHighValue>>8);
    sRandomRangeHighValue += sRandomRangeLowValue;
    sRandomRangeLowValue += sRandomRangeHighValue;
    int val = std::abs(sRandomRangeHighValue);
    return double(val) / double(RAND_MAX);
  }

  /// Returns the angle (given in radians) in degrees.
  static double deg(double a)  {
    return (180.0 / s_pi) * a;
  }

  /// Returns the angle (given in degrees) in radians.
  static double rad(double a)  {
    return  (s_pi / 180.0) * a;
  }

  /// Returns the minimum value among its arguments
  static double min(double x, double y) {
    return std::min(x, y);
  }

  /// Returns the maximum value among its arguments
  static double max(double x, double y) {
    return std::max(x, y);
  }

  /// Convert rectangular coordinates into polar radius.
  static double recttopolr(double x, double y) {
    return std::sqrt((x * x) + (y * y));
  }

  static double cosine_ramp(double t, double rampTime) {
    if( t < rampTime ) 
    {
      return (1.0 - std::cos(t*s_pi/rampTime))/2.0;
    }
    else 
    {
      return 1.0;
    }
  }

  /// Returns -1,0, or 1 depending on whether x is negative, zero, or positive.
  static double sign(double a)  {
    if(a==0.0) return 0.0;
    return (a > 0.0 ) ? 1.0 : -1.0;
  }

  /// Returns 1.0 if the input value is greater than the step time t.
  static double unit_step(double a, double t)  {
    return (a > t) ? 0.0 : -1.0;
  }

  /// Convert rectangular coordinates into polar angle.
  static double recttopola(double x, double y) {
    double tmp = std::atan2(y, x);

    /* Convert to 0.0 to 2 * PI */
    if (tmp < 0.0) {
      return tmp + (2.0 * s_pi);
    } else {
      return tmp;
    }
  }

  /// Convert polar coordinates (r,theta) into x coordinate.
  static double poltorectx(double r, double theta) {
    return r * std::cos(theta);
  }

  /// Convert polar coordinates (r,theta) into y coordinate.
  static double poltorecty(double r, double theta) {
    return r * std::sin(theta);
  }
}


CFunctionMap::CFunctionMap() 
{
  /// These random number functions support calls to
  /// the ANSI C random number generator.
  (*this)["rand"]        = new CFunction0(real_rand);
  (*this)["srand"]       = new CFunction1(real_srand);
  (*this)["randomize"]   = new CFunction0(randomize);

  /// These random number functions support a non-platform
  /// specific random number function.
  (*this)["random"]       = new CFunction0(random);
  (*this)["srandom"]      = new CFunction1(random_seed);

  (*this)["exp"]          = new CFunction1(std::exp);
  (*this)["ln"]           = new CFunction1(std::log);
  (*this)["log"]          = new CFunction1(std::log);
  (*this)["log10"]        = new CFunction1(std::log10);
  (*this)["pow"]          = new CFunction2(std::pow);
  (*this)["sqrt"]         = new CFunction1(std::sqrt);
  (*this)["erfc"]         = new CFunction1(erfc);
  (*this)["erf"]          = new CFunction1(erf);

  (*this)["acos"]         = new CFunction1(std::acos);
  (*this)["asin"]         = new CFunction1(std::asin);
  (*this)["atan"]         = new CFunction1(std::atan);
  (*this)["atan2"]        = new CFunction2(std::atan2);
  (*this)["ceil"]         = new CFunction1(std::ceil);
  (*this)["cos"]          = new CFunction1(std::cos);
  (*this)["cosh"]         = new CFunction1(std::cosh);
  (*this)["floor"]        = new CFunction1(std::floor);
  (*this)["sin"]          = new CFunction1(std::sin);
  (*this)["sinh"]         = new CFunction1(std::sinh);
  (*this)["tan"]          = new CFunction1(std::tan);
  (*this)["tanh"]         = new CFunction1(std::tanh);

  (*this)["abs"]          = new CFunction1(std::fabs);
  (*this)["fabs"]         = new CFunction1(std::fabs);
  (*this)["deg"]          = new CFunction1(deg);
  (*this)["mod"]          = new CFunction2(std::fmod);
  (*this)["fmod"]         = new CFunction2(std::fmod);
  (*this)["ipart"]        = new CFunction1(ipart);
  (*this)["fpart"]        = new CFunction1(fpart);
  (*this)["max"]          = new CFunction2(max);
  (*this)["min"]          = new CFunction2(min);
  (*this)["poltorectx"]   = new CFunction2(poltorectx);
  (*this)["poltorecty"]   = new CFunction2(poltorecty);
  (*this)["rad"]          = new CFunction1(rad);
  (*this)["recttopola"]   = new CFunction2(recttopola);
  (*this)["recttopolr"]   = new CFunction2(recttopolr);

  (*this)["cosine_ramp"]  = new CFunction2(cosine_ramp);
  (*this)["sign"]         = new CFunction1(sign);
  (*this)["unit_step"]    = new CFunction2(unit_step);
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
