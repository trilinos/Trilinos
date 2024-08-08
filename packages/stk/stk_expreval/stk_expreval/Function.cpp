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

#include "stk_expreval/Function.hpp"
#include <cmath>
#include <ctime>

namespace stk {
namespace expreval {

int sRandomRangeHighValue = 3191613;
int sRandomRangeLowValue  = 1739623;
typedef double (*CExtern0)();
typedef double (*CExtern1)(double);
typedef double (*CExtern2)(double, double);
typedef double (*CExtern3)(double, double, double);
typedef double (*CExtern4)(double, double, double, double);
typedef double (*CExtern5)(double, double, double, double, double);
typedef double (*CExtern8)(double, double, double, double, double, double, double, double);

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

template <>
class CFunction<CExtern5> : public CFunctionBase
{
public:
  typedef CExtern5 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(5),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv)
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 5 arguments"); }
#endif
    return (*m_function)(argv[0], argv[1], argv[2], argv[3], argv[4]);
  }

private:
  Signature m_function;
};

template <>
class CFunction<CExtern8> : public CFunctionBase
{
public:
  typedef CExtern8 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(8),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual double operator()(int argc, const double *argv)
  {
#ifndef NDEBUG
    if (argc != getArgCount()) { throw std::runtime_error("Argument count mismatch, function should have 8 arguments"); }
#endif
    return (*m_function)(argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
  }

private:
  Signature m_function;
};

typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;
typedef CFunction<CExtern3> CFunction3;
typedef CFunction<CExtern4> CFunction4;
typedef CFunction<CExtern5> CFunction5;
typedef CFunction<CExtern8> CFunction8;


CFunctionMap::CFunctionMap()
{
  (*this).emplace("rand",         new CFunction0(real_rand));
  (*this).emplace("srand",        new CFunction1(real_srand));

  (*this).emplace("random",          new CFunction0(random0));
  (*this).emplace("random",          new CFunction1(random1));
  (*this).emplace("time",            new CFunction0(current_time));
  (*this).emplace("ts_random",       new CFunction4(time_space_random));
  (*this).emplace("ts_normal",       new CFunction8(time_space_normal));

  (*this).emplace("exp",             new CFunction1(std::exp));
  (*this).emplace("ln",              new CFunction1(std::log));
  (*this).emplace("log",             new CFunction1(std::log));
  (*this).emplace("log10",           new CFunction1(std::log10));
  (*this).emplace("pow",             new CFunction2(std::pow));
  (*this).emplace("sqrt",            new CFunction1(std::sqrt));
  (*this).emplace("erfc",            new CFunction1(std::erfc));
  (*this).emplace("erf",             new CFunction1(std::erf));

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

  (*this).emplace("point2d",         new CFunction4(point_2));
  (*this).emplace("point3d",         new CFunction5(point_3));

  (*this).emplace("cos_ramp",        new CFunction1(cosine_ramp1));
  (*this).emplace("cos_ramp",        new CFunction2(cosine_ramp2));
  (*this).emplace("cos_ramp",        new CFunction3(cosine_ramp3));
  (*this).emplace("cosine_ramp",     new CFunction1(cosine_ramp1));
  (*this).emplace("cosine_ramp",     new CFunction2(cosine_ramp2));
  (*this).emplace("cosine_ramp",     new CFunction3(cosine_ramp3));
  (*this).emplace("linear_ramp",     new CFunction3(linear_ramp3));
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
