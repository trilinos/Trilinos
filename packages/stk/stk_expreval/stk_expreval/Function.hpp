// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef stk_expreval_function_hpp
#define stk_expreval_function_hpp

#include "Kokkos_Core.hpp"
#include "stk_expreval/Constants.hpp"
#include "stk_util/util/string_case_compare.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <string>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <ctime>
#include <iostream>

namespace stk {
namespace expreval {

enum class FunctionType {
  ABS,
  MAX,
  MIN,
  SIGN,
  IPART,
  FPART,
  CEIL,
  FLOOR,
  MOD,
  POW,
  SQRT,
  EXP,
  LN,
  LOG10,

  DEG,
  RAD,
  SIN,
  COS,
  TAN,
  ASIN,
  ACOS,
  ATAN,
  ATAN2,
  SINH,
  COSH,
  TANH,
  ASINH,
  ACOSH,
  ATANH,
  ERF,
  ERFC,
  POLTORECTX,
  POLTORECTY,
  RECTTOPOLR,
  RECTTOPOLA,

  UNIT_STEP,
  CYCLOIDAL_RAMP,
  COS_RAMP,
  HAVERSINE_PULSE,
  POINT2D,
  POINT3D,

  EXPONENTIAL_PDF,
  LOG_UNIFORM_PDF,
  NORMAL_PDF,
  WEIBULL_PDF,
  GAMMA_PDF,

  RAND,
  SRAND,
  RANDOM,
  TS_RANDOM,
  TS_NORMAL,
  TIME,

  UNDEFINED
};

KOKKOS_INLINE_FUNCTION
double cycloidal_ramp(double t, double t1, double t2)
{
  if (t < t1) {
    return 0.0;
  }
  else if (t < t2) {
    return (t-t1)/(t2-t1)-1/(two_pi())*sin(two_pi()/(t2-t1)*(t-t1));
  }
  else {
    return 1.0;
  }
}

/// extract signed integral value from floating-point number
KOKKOS_INLINE_FUNCTION
double ipart(double x)
{
  double y;
  std::modf(x, &y);
  return y;
}

/// Extract fractional value from floating-point number
KOKKOS_INLINE_FUNCTION
double fpart(double x)
{
  double y;
  return std::modf(x, &y);
}

/// Interface to the pseudo-random number generator function rand
/// provided by ANSI C math library.
KOKKOS_INLINE_FUNCTION
double real_rand()
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  return static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX) + 1.0);
#else
  NGP_ThrowErrorMsg("The rand function is not supported on GPUs");
  return 0.0;
#endif
}

/// Sets x as the random number seed. Interface to the srand function provided by the
/// ANSI C math library.
KOKKOS_INLINE_FUNCTION
double real_srand(double x)
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  std::srand(static_cast<int>(x));
  return 0.0;
#else
  NGP_ThrowErrorMsg("The srand function is not supported on GPUs");
  return 0.0;
#endif
}

/// Return the current time
KOKKOS_INLINE_FUNCTION
double current_time()
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  return static_cast<double>(::time(nullptr));
#else
  NGP_ThrowErrorMsg("The time function is not supported on GPUs");
  return 0.0;
#endif
}

extern int sRandomRangeHighValue;
extern int sRandomRangeLowValue;

/// Sets x as the "seed" for the pseudo-random number generator.
KOKKOS_INLINE_FUNCTION
void random_seed(double x)
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  int y = std::hash<double>{}(x);
  sRandomRangeHighValue =  y;
  sRandomRangeLowValue  = ~y;
#endif
}

/// Non-platform specific (pseudo) random number generator.
KOKKOS_INLINE_FUNCTION
double random0()
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  sRandomRangeHighValue = (sRandomRangeHighValue<<8) + (sRandomRangeHighValue>>8);
  sRandomRangeHighValue += sRandomRangeLowValue;
  sRandomRangeLowValue += sRandomRangeHighValue;
  int val = std::abs(sRandomRangeHighValue);
  return double(val) / double(RAND_MAX);
#else
  NGP_ThrowErrorMsg("The random function is not supported on GPUs");
  return 0.0;
#endif
}

/// Non-platform specific (pseudo) random number generator.
KOKKOS_INLINE_FUNCTION
double random1(double seed)
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  random_seed(seed);
  return random0();
#else
  NGP_ThrowErrorMsg("The random function is not supported on GPUs");
  return 0.0;
#endif
}

/// Non-platform specific (pseudo) random number generator that
/// is deterministic for a given point in time and space
KOKKOS_INLINE_FUNCTION
double time_space_random(double t, double x, double y, double z)
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  double ts = t + x + y + z + x*y + y*z + x*z + x*y*z;
  random_seed(ts);
  return random0();
#else
  NGP_ThrowErrorMsg("The ts_random function is not supported on GPUs");
  return 0.0;
#endif
}

KOKKOS_INLINE_FUNCTION
double time_space_normal(double t, double x, double y, double z, double mu, double sigma, double minR, double maxR)
{
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  double ts = t + x + y + z + x*y + y*z + x*z + x*y*z;
  random_seed(ts);

  static const double epsilon = DBL_MIN;

  // Box-Muller transformation from two uniform random numbers
  // to a gaussian distribution
  double u1 = std::fmax(epsilon, random0());
  double u2 = std::fmax(epsilon, random0());

  double z0 = std::sqrt(-2.0 * std::log(u1)) * std::cos(two_pi() * u2);

  return std::fmax(minR, std::fmin(maxR, z0*sigma + mu));
#else
  NGP_ThrowErrorMsg("The ts_normal function is not supported on GPUs");
  return 0.0;
#endif
}

/// Returns the angle (input in radians) in degrees.
KOKKOS_INLINE_FUNCTION
double deg(double a)
{
  return radian_to_degree() * a;
}

/// Returns the angle (input in degrees) in radians.
KOKKOS_INLINE_FUNCTION
double rad(double a)
{
  return degree_to_radian() * a;
}

/// Returns the minimum value among its arguments
KOKKOS_INLINE_FUNCTION
double min_2(double a, double b)
{
  return std::fmin(a, b);
}

/// Returns the minimum value among its arguments
KOKKOS_INLINE_FUNCTION
double min_3(double a, double b, double c)
{
  return std::fmin(std::fmin(a, b), c);
}

/// Returns the minimum value among its arguments
KOKKOS_INLINE_FUNCTION
double min_4(double a, double b, double c, double d)
{
  return std::fmin(std::fmin(a, b), std::fmin(c,d));
}

/// Returns the maximum value among its arguments
KOKKOS_INLINE_FUNCTION
double max_2(double a, double b)
{
  return std::fmax(a, b);
}

/// Returns the maximum value among its arguments
KOKKOS_INLINE_FUNCTION
double max_3(double a, double b, double c)
{
  return std::fmax(std::fmax(a, b), c);
}

/// Returns the maximum value among its arguments
KOKKOS_INLINE_FUNCTION
double max_4(double a, double b, double c, double d)
{
  return std::fmax(std::fmax(a, b), std::fmax(c,d));
}

/// Convert rectangular coordinates into polar radius.
KOKKOS_INLINE_FUNCTION
double recttopolr(double x, double y)
{
  return std::sqrt((x * x) + (y * y));
}

KOKKOS_INLINE_FUNCTION
double cosine_ramp3(double t, double t1, double t2)
{
  if (t < t1) {
    return 0.0;
  }
  else if (t < t2) {
    return (1.0 - std::cos((t-t1)*pi() /(t2-t1)))/2.0;
  }
  else {
    return 1.0;
  }
}

KOKKOS_INLINE_FUNCTION
double haversine_pulse(double t, double t1, double t2)
{
  if (t < t1) {
    return 0.0;
  }
  else if (t < t2) {
    return std::pow(std::sin(pi() *(t-t1)/(t2-t1)),2);
  }
  else {
    return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
double point_2(double x, double y, double r, double w)
{
  const double ri = std::sqrt(x*x + y*y);
  return 1.0 - cosine_ramp3(ri, r-0.5*w, r+0.5*w);
}

KOKKOS_INLINE_FUNCTION
double point_3(double x, double y, double z, double r, double w)
{
  const double ri = std::sqrt(x*x + y*y + z*z);
  return 1.0 - cosine_ramp3(ri, r-0.5*w, r+0.5*w);
}

KOKKOS_INLINE_FUNCTION
double cosine_ramp1(double t)
{
  return cosine_ramp3(t, 0.0, 1.0);
}

KOKKOS_INLINE_FUNCTION
double cosine_ramp2(double t, double rampEndTime)
{
  return cosine_ramp3(t, 0.0, rampEndTime);
}

/// Weibull distribution probability distribution function.
KOKKOS_INLINE_FUNCTION
double weibull_pdf(double x, double shape, double scale)
{
  return (x >= 0) ? (shape/scale)*std::pow(x/scale, shape-1)*std::exp(-std::pow(x/scale, shape)) : 0;
}

/// Normal (Gaussian) distribution probability distribution function.
KOKKOS_INLINE_FUNCTION
double normal_pdf(double x, double mean, double standard_deviation)
{
  return std::exp(-(x-mean)*(x-mean)/(2.0*standard_deviation*standard_deviation)) /
         std::sqrt(2.0*pi()*standard_deviation*standard_deviation);
}

/// Exponential Uniform distribution probability distribution function
KOKKOS_INLINE_FUNCTION
double exponential_pdf(double x, double beta)
{
  return (x >= 0.0) ? std::exp(-x/beta)/beta : 0.0;
}

/// Log Uniform distribution probability distribution function
KOKKOS_INLINE_FUNCTION
double log_uniform_pdf(double x, double lower_range, double upper_range)
{
  return (x >= lower_range && x <= upper_range) ? 1.0/((std::log(upper_range) - std::log(lower_range))*x) : 0.0;
}

/// Gamma continuous probability distribution function.
KOKKOS_INLINE_FUNCTION
double gamma_pdf(double x, double shape, double scale)
{
  return (x >= 0) ? 1/(std::tgamma(shape)*std::pow(scale, shape))*std::pow(x, shape-1)*std::exp(-x/scale) : 0;
}

/// Returns -1 or 1 depending on whether x is negative or positive.
KOKKOS_INLINE_FUNCTION
double sign(double a)
{
  return (a >= 0.0) ? 1.0 : -1.0;
}

/// Returns 1.0 if the input value t is greater than tstart and less than tstop.
KOKKOS_INLINE_FUNCTION
double unit_step3(double t, double tstart, double tstop)
{
  return (t < tstart || t > tstop) ? 0.0 : 1.0;
}

/// Convert rectangular coordinates into polar angle.
KOKKOS_INLINE_FUNCTION
double recttopola(double x, double y)
{
  double tmp = std::atan2(y, x);
  // Convert to 0.0 to 2 * PI
  return (tmp < 0.0) ? tmp + two_pi() : tmp;
}

/// Convert polar coordinates (r,theta) into x coordinate.
KOKKOS_INLINE_FUNCTION
double poltorectx(double r, double theta)
{
  return r * std::cos(theta);
}

/// Convert polar coordinates (r,theta) into y coordinate.
KOKKOS_INLINE_FUNCTION
double poltorecty(double r, double theta)
{
  return r * std::sin(theta);
}

class CFunctionBase
{
public:
  explicit CFunctionBase(int arg_count)
    : m_argCount(arg_count)
  {}

  virtual ~CFunctionBase()
  {}

  virtual double operator()(int argc, const double * argv) = 0;

  int getArgCount() const { return m_argCount; }

private:
  int  m_argCount;
};


template <class S>
class CFunction;

class CFunctionMap : public std::multimap<std::string, CFunctionBase *, LessCase>
{
public:
  CFunctionMap();
  ~CFunctionMap();
};

CFunctionMap &getCFunctionMap();

inline void addFunction(const std::string & name, CFunctionBase * function) {
  getCFunctionMap().insert(std::make_pair(name, function));
}

} // namespace expreval
} // namespace stk

#endif // stk_expreval_function_hpp
