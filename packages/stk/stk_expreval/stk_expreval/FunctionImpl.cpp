#include <math.h>
#include <cmath>
#include <ctime>
#include <math.h>       //Needed for erf and erfc on solaris.

#include <stk_expreval/Function.hpp>
#include <stk_expreval/Constants.hpp>

#include <boost/math/distributions.hpp>

namespace stk {
namespace expreval {

  namespace bmp  = boost::math::policies;

typedef boost::math::
  weibull_distribution< double,
                       boost::math::policies::policy< bmp::overflow_error<bmp::ignore_error> > >
  weibull_dist;

typedef boost::math::
  gamma_distribution< double,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  gamma_dist;

typedef boost::math::
  normal_distribution< double,
                       bmp::policy< bmp::overflow_error<bmp::ignore_error> > >
  normal_dist;

extern "C" {
  typedef double (*CExtern0)();
  typedef double (*CExtern1)(double);
  typedef double (*CExtern2)(double, double);
  typedef double (*CExtern3)(double, double, double);
  typedef double (*CExtern4)(double, double, double, double);
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
      throw std::runtime_error("Argument count mismatch, function should have 0 arguments");

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
      throw std::runtime_error("Argument count mismatch, function should have 1 argument");

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
      throw std::runtime_error("Argument count mismatch, function should have 2 arguments");

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
      throw std::runtime_error("Argument count mismatch, function should have 3 arguments");

    return (*m_function)(argv[0], argv[1], argv[2]);
  }

private:
  Signature	m_function;
};

template <>
class CFunction<CExtern4> : public CFunctionBase
{
public:
  typedef CExtern4 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(4),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv) {
    // KHP: Maybe only check in debug?
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch, function should have 4 arguments");
    return (*m_function)(argv[0], argv[1], argv[2], argv[3]);
  }

private:
  Signature	m_function;
};

typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;
typedef CFunction<CExtern3> CFunction3;
typedef CFunction<CExtern4> CFunction4;


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

  /// Return the current time
  static double current_time() {
    return static_cast<double>(::time(NULL));
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
  static double random0() {
    sRandomRangeHighValue = (sRandomRangeHighValue<<8) + (sRandomRangeHighValue>>8);
    sRandomRangeHighValue += sRandomRangeLowValue;
    sRandomRangeLowValue += sRandomRangeHighValue;
    int val = std::abs(sRandomRangeHighValue);
    return double(val) / double(RAND_MAX);
  }

  /// Non-platform specific (pseudo) random number generator.
  static double random1(double seed) {
    random_seed(seed);
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
  static double min2(double a, double b) {
    return std::min(a, b);
  }

  /// Returns the minimum value among its arguments
  static double min3(double a, double b, double c) {
    return std::min(std::min(a, b), c);
  }

  /// Returns the minimum value among its arguments
  static double min4(double a, double b, double c, double d) {
    return std::min(std::min(a, b), std::min(c,d));
  }

  /// Returns the maximum value among its arguments
  static double max2(double a, double b) {
    return std::max(a, b);
  }

  /// Returns the maximum value among its arguments
  static double max3(double a, double b, double c) {
    return std::max(std::max(a, b), c);
  }

  /// Returns the maximum value among its arguments
  static double max4(double a, double b, double c, double d) {
    return std::max(std::max(a, b), std::max(c,d));
  }

  /// Convert rectangular coordinates into polar radius.
  static double recttopolr(double x, double y) {
    return std::sqrt((x * x) + (y * y));
  }

  static double cosine_ramp3(double t, double rampStartTime, double rampEndTime) {
    if( t < rampStartTime    )
    {
      return 0.0;
    }
    else if( t < rampEndTime )
    {
      return (1.0 - std::cos((t-rampStartTime)*s_pi/(rampEndTime-rampStartTime)))/2.0;
    }
    else 
    {
      return 1.0;
    }
  }

  static double cosine_ramp1(double t) {
    return cosine_ramp3(t, 0.0, 1.0);
  }

  static double cosine_ramp2(double t, double rampEndTime) {
    return cosine_ramp3(t, 0.0, rampEndTime);
  }

  /// Weibull distribution probability distribution function.
  double weibull_pdf(double x, double shape, double scale)
  {
#if defined(__PATHSCALE__)
    return 0.0;
#else
    weibull_dist weibull1(shape, scale);
    return boost::math::pdf(weibull1, x);
#endif
  }

  /// Normal (Gaussian) distribution probability distribution function.
  double normal_pdf(double x, double mean, double standard_deviation)
  {
#if defined(__PATHSCALE__)
    return 0.0;
#else
    normal_dist normal1(mean, standard_deviation);
    return boost::math::pdf(normal1, x);
#endif
  }

  /// Uniform distribution probability distribution function.
  double uniform_pdf(double lower_range, double upper_range)
  {
    // Note, no error checking here...
    return 1.0/(upper_range - lower_range);
  }
  
  /// Exponential Uniform distribution probability distribution function
  inline double exponential_pdf(double x, double beta)
  { return std::exp(-x/beta)/beta; }

  /// Log Uniform distribution probability distribution function
  inline double log_uniform_pdf(double x, double lower_range, double upper_range)
  { return 1.0/(std::log(upper_range) - std::log(lower_range))/x; }

  /// Gamma continuous probability distribution function.
  inline double gamma_pdf(double x, double shape, double scale)
  {
#if defined(__PATHSCALE__)
    return 0.0;
#else
    return boost::math::pdf(gamma_dist(shape,scale), x);
#endif
  }

  inline double phi(double beta)
  {
#if defined(__PATHSCALE__)
    return 0.0;
#else
    return boost::math::pdf(normal_dist(0.,1.), beta);
#endif
  }

  /// Returns a probability < 0.5 for negative beta and a probability > 0.5 for positive beta.
  inline double Phi(double beta)
  {
#if defined(__PATHSCALE__)
    return 0.0;
#else
    return boost::math::cdf(normal_dist(0.,1.), beta);
#endif
  }

  inline double bounded_normal_pdf(double x, double mean, double std_dev, double lwr, double upr)
  {
    double Phi_lms = (lwr > -std::numeric_limits<double>::max()) ? Phi((lwr-mean)/std_dev) : 0.;
    double Phi_ums = (upr <  std::numeric_limits<double>::max()) ? Phi((upr-mean)/std_dev) : 1.;
    return phi((x-mean)/std_dev)/(Phi_ums - Phi_lms)/std_dev;
  }

  /// Returns -1 or 1 depending on whether x is negative or positive.
  static double sign(double a)  {
    return (a >= 0.0 ) ? 1.0 : -1.0;
  }

  /// Returns 1.0 if the input value t is greater than tstart and less than tstop.
  static double unit_step3(double t, double tstart, double tstop)  {
    return (t < tstart || t > tstop) ? 0.0 : 1.0;
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
  (*this).insert(std::make_pair("rand",         new CFunction0(real_rand)));
  (*this).insert(std::make_pair("srand",        new CFunction1(real_srand)));
  (*this).insert(std::make_pair("randomize",    new CFunction0(randomize)));

  /// These random number functions support a non-platform
  /// specific random number function.
  (*this).insert(std::make_pair("time",            new CFunction0(current_time)));
  (*this).insert(std::make_pair("random",          new CFunction0(random0)));
  (*this).insert(std::make_pair("random",          new CFunction1(random1)));

  (*this).insert(std::make_pair("exp",             new CFunction1(std::exp)));
  (*this).insert(std::make_pair("ln",              new CFunction1(std::log)));
  (*this).insert(std::make_pair("log",             new CFunction1(std::log)));
  (*this).insert(std::make_pair("log10",           new CFunction1(std::log10)));
  (*this).insert(std::make_pair("pow",             new CFunction2(std::pow)));
  (*this).insert(std::make_pair("sqrt",            new CFunction1(std::sqrt)));
  (*this).insert(std::make_pair("erfc",            new CFunction1(erfc)));
  (*this).insert(std::make_pair("erf",             new CFunction1(erf)));

  (*this).insert(std::make_pair("acos",            new CFunction1(std::acos)));
  (*this).insert(std::make_pair("asin",            new CFunction1(std::asin)));
  (*this).insert(std::make_pair("atan",            new CFunction1(std::atan)));
  (*this).insert(std::make_pair("atan2",           new CFunction2(std::atan2)));
  (*this).insert(std::make_pair("ceil",            new CFunction1(std::ceil)));
  (*this).insert(std::make_pair("cos",             new CFunction1(std::cos)));
  (*this).insert(std::make_pair("cosh",            new CFunction1(std::cosh)));
  (*this).insert(std::make_pair("floor",           new CFunction1(std::floor)));
  (*this).insert(std::make_pair("sin",             new CFunction1(std::sin)));
  (*this).insert(std::make_pair("sinh",            new CFunction1(std::sinh)));
  (*this).insert(std::make_pair("tan",             new CFunction1(std::tan)));
  (*this).insert(std::make_pair("tanh",            new CFunction1(std::tanh)));

  (*this).insert(std::make_pair("abs",             new CFunction1(std::fabs)));
  (*this).insert(std::make_pair("fabs",            new CFunction1(std::fabs)));
  (*this).insert(std::make_pair("deg",             new CFunction1(deg)));
  (*this).insert(std::make_pair("mod",             new CFunction2(std::fmod)));
  (*this).insert(std::make_pair("fmod",            new CFunction2(std::fmod)));
  (*this).insert(std::make_pair("ipart",           new CFunction1(ipart)));
  (*this).insert(std::make_pair("fpart",           new CFunction1(fpart)));
  (*this).insert(std::make_pair("max",             new CFunction2(max2)));
  (*this).insert(std::make_pair("max",             new CFunction3(max3)));
  (*this).insert(std::make_pair("max",             new CFunction4(max4)));
  (*this).insert(std::make_pair("min",             new CFunction2(min2)));
  (*this).insert(std::make_pair("min",             new CFunction3(min3)));
  (*this).insert(std::make_pair("min",             new CFunction4(min4)));
  (*this).insert(std::make_pair("poltorectx",      new CFunction2(poltorectx)));
  (*this).insert(std::make_pair("poltorecty",      new CFunction2(poltorecty)));
  (*this).insert(std::make_pair("rad",             new CFunction1(rad)));
  (*this).insert(std::make_pair("recttopola",      new CFunction2(recttopola)));
  (*this).insert(std::make_pair("recttopolr",      new CFunction2(recttopolr)));

  (*this).insert(std::make_pair("cos_ramp",        new CFunction1(cosine_ramp1)));
  (*this).insert(std::make_pair("cos_ramp",        new CFunction2(cosine_ramp2)));
  (*this).insert(std::make_pair("cos_ramp",        new CFunction3(cosine_ramp3)));
  (*this).insert(std::make_pair("cosine_ramp",     new CFunction1(cosine_ramp1)));
  (*this).insert(std::make_pair("cosine_ramp",     new CFunction2(cosine_ramp2)));
  (*this).insert(std::make_pair("cosine_ramp",     new CFunction3(cosine_ramp3)));

  (*this).insert(std::make_pair("sign",            new CFunction1(sign)));
  (*this).insert(std::make_pair("unit_step",       new CFunction3(unit_step3)));

  (*this).insert(std::make_pair("weibull_pdf",     new CFunction3(weibull_pdf)));
  (*this).insert(std::make_pair("normal_pdf",      new CFunction3(normal_pdf)));
  (*this).insert(std::make_pair("gamma_pdf",       new CFunction3(gamma_pdf)));
  (*this).insert(std::make_pair("log_uniform_pdf", new CFunction3(log_uniform_pdf)));
  (*this).insert(std::make_pair("uniform_pdf",     new CFunction2(uniform_pdf)));
  (*this).insert(std::make_pair("exponential_pdf", new CFunction2(exponential_pdf)));
}

CFunctionMap::~CFunctionMap()
{
  for (CFunctionMap::iterator it = begin(); it != end(); ++it) {
    delete (*it).second;
  }
}

CFunctionMap &
getCFunctionMap()
{
  static CFunctionMap s_functionMap;

  return s_functionMap;
}

} // namespace expreval
} // namespace stk
