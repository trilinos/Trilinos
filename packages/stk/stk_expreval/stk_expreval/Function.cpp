// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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

  using weibull_dist = boost::math::weibull_distribution< double, boost::math::policies::policy< bmp::overflow_error<bmp::ignore_error> > >;
  using gamma_dist   = boost::math::gamma_distribution< double, bmp::policy< bmp::overflow_error<bmp::ignore_error> > >;
  using normal_dist  = boost::math::normal_distribution< double, bmp::policy< bmp::overflow_error<bmp::ignore_error> > >;

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

  virtual double operator()(int argc, const double *argv) 
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 0 arguments"); }
#endif
    return (*m_function)();
  }

private:
  Signature m_function;
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

  virtual double operator()(int argc, const double *argv) 
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 1 argument"); }
#endif
    return (*m_function)(argv[0]);
  }

private:
  Signature m_function;
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

  virtual double operator()(int argc, const double *argv) 
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 2 arguments"); }
#endif
    return (*m_function)(argv[0], argv[1]);
  }

private:
  Signature m_function;
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

  virtual double operator()(int argc, const double *argv) 
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 3 arguments"); }
#endif
    return (*m_function)(argv[0], argv[1], argv[2]);
  }

private:
  Signature m_function;
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

  virtual double operator()(int argc, const double *argv) 
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 4 arguments"); }
#endif
    return (*m_function)(argv[0], argv[1], argv[2], argv[3]);
  }

private:
  Signature m_function;
};

typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;
typedef CFunction<CExtern3> CFunction3;
typedef CFunction<CExtern4> CFunction4;


extern "C" {
  double cycloidal_ramp(double t, double t1, double t2)
  {
    if( t < t1 )
    {
      return 0.0;
    }
    else if( t < t2 )
    {
      return (t-t1)/(t2-t1)-1/(two_pi())*sin(two_pi()/(t2-t1)*(t-t1));
    }
    else 
    {
      return 1.0;
    }
  }
}

namespace {
extern "C" {
  /// extract signed integral value from floating-point number
  double ipart(double x) 
  {
    double y;
    std::modf(x, &y);
    return y;
  }

  /// Extract fractional value from floating-point number
  double fpart(double x) 
  {
    double y;
    return std::modf(x, &y);
  }

  /// Interface to the pseudo-random number generator function rand
  /// provided by ANSI C math library.
  double real_rand() 
  {
    return static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX) + 1.0);
  }

  /// Sets x as the random number seed. Interface to the srand function provided by the
  /// ANSI C math library.
  double real_srand(double x) 
  {
    std::srand(static_cast<int>(x));
    return 0.0;
  }

  /// Return the current time
  double current_time()
  {
    return static_cast<double>(::time(nullptr));
  }

  /// Sets x as the "seed" for the pseudo-random number generator.
  void random_seed(double x) 
  {
    int y = static_cast<int>(x);
    sRandomRangeHighValue =  y;
    sRandomRangeLowValue  = ~y;
  }

  /// Non-platform specific (pseudo) random number generator.
  double random0() 
  {
    sRandomRangeHighValue = (sRandomRangeHighValue<<8) + (sRandomRangeHighValue>>8);
    sRandomRangeHighValue += sRandomRangeLowValue;
    sRandomRangeLowValue += sRandomRangeHighValue;
    int val = std::abs(sRandomRangeHighValue);
    return double(val) / double(RAND_MAX);
  }

  /// Non-platform specific (pseudo) random number generator.
  double random1(double seed) 
  {
    random_seed(seed);
    return random0();
  }

  /// Returns the angle (input in radians) in degrees.
  double deg(double a)  
  {
    return radian_to_degree() * a;
  }

  /// Returns the angle (input in degrees) in radians.
  double rad(double a)  
  {
    return degree_to_radian() * a;
  }

  /// Returns the minimum value among its arguments
  double min_2(double a, double b) 
  {
    return std::min(a, b);
  }

  /// Returns the minimum value among its arguments
  double min_3(double a, double b, double c) 
  {
    return std::min(std::min(a, b), c);
  }

  /// Returns the minimum value among its arguments
  double min_4(double a, double b, double c, double d) 
  {
    return std::min(std::min(a, b), std::min(c,d));
  }

  /// Returns the maximum value among its arguments
  double max_2(double a, double b) 
  {
    return std::max(a, b);
  }

  /// Returns the maximum value among its arguments
  double max_3(double a, double b, double c) 
  {
    return std::max(std::max(a, b), c);
  }

  /// Returns the maximum value among its arguments
  double max_4(double a, double b, double c, double d) 
  {
    return std::max(std::max(a, b), std::max(c,d));
  }

  /// Convert rectangular coordinates into polar radius.
  double recttopolr(double x, double y) 
  {
    return std::sqrt((x * x) + (y * y));
  }

  double cosine_ramp3(double t, double t1, double t2) 
  {
    if( t < t1    )
    {
      return 0.0;
    }
    else if( t < t2 )
    {
      return (1.0 - std::cos((t-t1)*pi() /(t2-t1)))/2.0;
    }
    else 
    {
      return 1.0;
    }
  }

  double haversine_pulse(double t, double t1, double t2)
  {
    if( t < t1 )
    {
      return 0.0;
    }
    else if( t < t2 )
    {
      return std::pow(std::sin(pi() *(t-t1)/(t2-t1)),2);
    }
    else 
    {
      return 0.0;
    }
  }

  double cosine_ramp1(double t) 
  {
    return cosine_ramp3(t, 0.0, 1.0);
  }

  double cosine_ramp2(double t, double rampEndTime) 
  {
    return cosine_ramp3(t, 0.0, rampEndTime);
  }

  /// Weibull distribution probability distribution function.
  double weibull_pdf(double x, double shape, double scale)
  {
    weibull_dist weibull1(shape, scale);
    return boost::math::pdf(weibull1, x);
  }

  /// Normal (Gaussian) distribution probability distribution function.
  double normal_pdf(double x, double mean, double standard_deviation)
  {
    normal_dist normal1(mean, standard_deviation);
    return boost::math::pdf(normal1, x);
  }

  /// Exponential Uniform distribution probability distribution function
  double exponential_pdf(double x, double beta)
  { 
    return std::exp(-x/beta)/beta; 
  }

  /// Log Uniform distribution probability distribution function
  double log_uniform_pdf(double x, double lower_range, double upper_range) 
  { 
    return 1.0/(std::log(upper_range) - std::log(lower_range))/x; 
  }

  /// Gamma continuous probability distribution function.
  double gamma_pdf(double x, double shape, double scale)
  {
    return boost::math::pdf(gamma_dist(shape,scale), x);
  }

  /// Returns -1 or 1 depending on whether x is negative or positive.
  double sign(double a)  
  {
    return ( a >= 0.0 ) ? 1.0 : -1.0;
  }

  /// Returns 1.0 if the input value t is greater than tstart and less than tstop.
  double unit_step3(double t, double tstart, double tstop)  
  {
    return (t < tstart || t > tstop) ? 0.0 : 1.0;
  }

  /// Convert rectangular coordinates into polar angle.
  double recttopola(double x, double y) 
  {
    double tmp = std::atan2(y, x);
    // Convert to 0.0 to 2 * PI
    return ( tmp < 0.0 ) ? tmp + two_pi() : tmp;
  }

  /// Convert polar coordinates (r,theta) into x coordinate.
  double poltorectx(double r, double theta) 
  {
    return r * std::cos(theta);
  }

  /// Convert polar coordinates (r,theta) into y coordinate.
  double poltorecty(double r, double theta) 
  {
    return r * std::sin(theta);
  }
}
}

CFunctionMap::CFunctionMap() 
{
  /// These random number functions support calls to
  /// the ANSI C random number generator.
  (*this).emplace("rand",         new CFunction0(real_rand));
  (*this).emplace("srand",        new CFunction1(real_srand));

  /// These random number functions support a platform
  /// independent random number function.
  (*this).emplace("random",          new CFunction0(random0));
  (*this).emplace("random",          new CFunction1(random1));
  (*this).emplace("time",            new CFunction0(current_time));

  (*this).emplace("exp",             new CFunction1(std::exp));
  (*this).emplace("ln",              new CFunction1(std::log));
  (*this).emplace("log",             new CFunction1(std::log));
  (*this).emplace("log10",           new CFunction1(std::log10));
  (*this).emplace("pow",             new CFunction2(std::pow));
  (*this).emplace("sqrt",            new CFunction1(std::sqrt));
  (*this).emplace("erfc",            new CFunction1(erfc));
  (*this).emplace("erf",             new CFunction1(erf));

  (*this).emplace("acos",            new CFunction1(std::acos));
  (*this).emplace("asin",            new CFunction1(std::asin));
  (*this).emplace("asinh",           new CFunction1(std::asinh));
  (*this).emplace("atan",            new CFunction1(std::atan));
  (*this).emplace("atan2",           new CFunction2(std::atan2));
  (*this).emplace("atanh",           new CFunction1(std::atanh));
  (*this).emplace("ceil",            new CFunction1(std::ceil));
  (*this).emplace("cos",             new CFunction1(std::cos));
  (*this).emplace("cosh",            new CFunction1(std::cosh));
  (*this).emplace("acosh",           new CFunction1(std::acosh));
  (*this).emplace("floor",           new CFunction1(std::floor));
  (*this).emplace("sin",             new CFunction1(std::sin));
  (*this).emplace("sinh",            new CFunction1(std::sinh));
  (*this).emplace("tan",             new CFunction1(std::tan));
  (*this).emplace("tanh",            new CFunction1(std::tanh));

  (*this).emplace("abs",             new CFunction1(std::fabs));
  (*this).emplace("fabs",            new CFunction1(std::fabs));
  (*this).emplace("deg",             new CFunction1(deg));
  (*this).emplace("mod",             new CFunction2(std::fmod));
  (*this).emplace("fmod",            new CFunction2(std::fmod));
  (*this).emplace("ipart",           new CFunction1(ipart));
  (*this).emplace("fpart",           new CFunction1(fpart));
  (*this).emplace("max",             new CFunction2(max_2));
  (*this).emplace("max",             new CFunction3(max_3));
  (*this).emplace("max",             new CFunction4(max_4));
  (*this).emplace("min",             new CFunction2(min_2));
  (*this).emplace("min",             new CFunction3(min_3));
  (*this).emplace("min",             new CFunction4(min_4));
  (*this).emplace("poltorectx",      new CFunction2(poltorectx));
  (*this).emplace("poltorecty",      new CFunction2(poltorecty));
  (*this).emplace("rad",             new CFunction1(rad));
  (*this).emplace("recttopola",      new CFunction2(recttopola));
  (*this).emplace("recttopolr",      new CFunction2(recttopolr));

  (*this).emplace("cos_ramp",        new CFunction1(cosine_ramp1));
  (*this).emplace("cos_ramp",        new CFunction2(cosine_ramp2));
  (*this).emplace("cos_ramp",        new CFunction3(cosine_ramp3));
  (*this).emplace("cosine_ramp",     new CFunction1(cosine_ramp1));
  (*this).emplace("cosine_ramp",     new CFunction2(cosine_ramp2));
  (*this).emplace("cosine_ramp",     new CFunction3(cosine_ramp3));
  (*this).emplace("haversine_pulse", new CFunction3(haversine_pulse));
  (*this).emplace("cycloidal_ramp",  new CFunction3(cycloidal_ramp));

  (*this).emplace("sign",            new CFunction1(sign));
  (*this).emplace("unit_step",       new CFunction3(unit_step3));

  (*this).emplace("weibull_pdf",     new CFunction3(weibull_pdf));
  (*this).emplace("normal_pdf",      new CFunction3(normal_pdf));
  (*this).emplace("gamma_pdf",       new CFunction3(gamma_pdf));
  (*this).emplace("log_uniform_pdf", new CFunction3(log_uniform_pdf));
  (*this).emplace("exponential_pdf", new CFunction2(exponential_pdf));
}

CFunctionMap::~CFunctionMap()
{
  for (CFunctionMap::iterator it = begin(); it != end(); ++it) 
  {
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
