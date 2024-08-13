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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_expreval/Evaluator.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <memory>

namespace {
//-BEGIN
double evaluate(const std::string & expression)
{
  stk::expreval::Eval eval(expression);
  eval.parse();
  return eval.evaluate();
}

TEST(HostEvaluation, testOpcodes)
{
  EXPECT_DOUBLE_EQ(evaluate(""), 0.0);
  EXPECT_DOUBLE_EQ(evaluate("-1.333"), -1.333);
  EXPECT_DOUBLE_EQ(evaluate("(1+4+9+16)+(9+16+25+36)"), 116);
  EXPECT_DOUBLE_EQ(evaluate("(-1-4)-(9-16)"), 2);
  EXPECT_DOUBLE_EQ(evaluate("2*3*4*5"),      120);
  EXPECT_DOUBLE_EQ(evaluate("(120/5)/(4/3)"), 18);
  EXPECT_DOUBLE_EQ(evaluate("---2"),         -2);
  EXPECT_DOUBLE_EQ(evaluate("9.1 % 4"),     1.1);
  EXPECT_DOUBLE_EQ(evaluate("3^2^2"),        81);
  EXPECT_DOUBLE_EQ(evaluate("0.1==0.999999"), 0);
  EXPECT_DOUBLE_EQ(evaluate("2!=(1+1)"),      0);
  EXPECT_DOUBLE_EQ(evaluate("1<1.000001"), 1);
  EXPECT_DOUBLE_EQ(evaluate("1>1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("1<=1"),       1);
  EXPECT_DOUBLE_EQ(evaluate("2>=(1+2)"),   0);
  EXPECT_ANY_THROW(evaluate("1 <= 2 < 3"));
  EXPECT_DOUBLE_EQ(evaluate("!0"),        1);
  EXPECT_DOUBLE_EQ(evaluate("!0.000001"), 0);
  EXPECT_DOUBLE_EQ(evaluate("0 && 0"),        0);
  EXPECT_DOUBLE_EQ(evaluate("0 && 1"),        0);
  EXPECT_DOUBLE_EQ(evaluate("0 || 0"),        0);
  EXPECT_DOUBLE_EQ(evaluate("0 || 1"),        1);
  EXPECT_DOUBLE_EQ(evaluate("0 ? 1 : (1+1)"),   2);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; z=3"),   3);
  EXPECT_DOUBLE_EQ(evaluate("x=1; y=2; x+y"),   3);
  EXPECT_DOUBLE_EQ(evaluate("(1+2+3+4)^(1+1)"), 100);
  EXPECT_DOUBLE_EQ(evaluate("15%(1+1+1)"),      0);
  EXPECT_DOUBLE_EQ(evaluate("x + y + z"),        0);
  EXPECT_DOUBLE_EQ(evaluate("x[0]"),             0);
  EXPECT_ANY_THROW(evaluate("x[0]+x[1]+x[2]"));
}

TEST(HostEvaluation, testFunctions)
{
  EXPECT_DOUBLE_EQ(evaluate("abs(-2*3)"), 6);
  EXPECT_DOUBLE_EQ(evaluate("fabs(1.5)"),  1.5);
  EXPECT_DOUBLE_EQ(evaluate("max(-1,-2,-3)"),    -1);
  EXPECT_DOUBLE_EQ(evaluate("min(3+2,2+1)"),      3);
  EXPECT_DOUBLE_EQ(evaluate("sign(-0.5)"),      -1);
  EXPECT_DOUBLE_EQ(evaluate("ipart(2.5)"),        2);
  EXPECT_DOUBLE_EQ(evaluate("fpart(-2.5)"),      -0.5);
  EXPECT_DOUBLE_EQ(evaluate("ceil(-0.999999)"),  0);
  EXPECT_DOUBLE_EQ(evaluate("floor(-0.000001)"), -1);
  EXPECT_DOUBLE_EQ(evaluate("mod(9, -4)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("fmod(9, 3.5)"),  2);
  EXPECT_DOUBLE_EQ(evaluate("pow(3, 2.5)"),               std::pow(3, 2.5));
  EXPECT_DOUBLE_EQ(evaluate("sqrt(1.21)"), 1.1);
  EXPECT_DOUBLE_EQ(evaluate("exp(1.5)"), std::exp(1.5));
  EXPECT_DOUBLE_EQ(evaluate("ln(exp(1.5))"), 1.5);
  EXPECT_DOUBLE_EQ(evaluate("log(0.5)"),      std::log(0.5));
  EXPECT_DOUBLE_EQ(evaluate("log10(1)"),        0);
  EXPECT_DOUBLE_EQ(evaluate("deg(PI/4)"),     45);
  EXPECT_DOUBLE_EQ(evaluate("rad(45)"),    stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("sin(PI/4)"),    std::sqrt(2)/2);
  EXPECT_DOUBLE_EQ(evaluate("cos(PI/4)"),    std::sqrt(2)/2);
  EXPECT_DOUBLE_EQ(evaluate("tan(PI/4)"),    1);
  EXPECT_DOUBLE_EQ(evaluate("asin(sqrt(2)/2)"), stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("acos(sqrt(2)/2)"), stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("atan(1)"),       stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("atan2(1, 1)"),    stk::expreval::pi()/4);
  EXPECT_DOUBLE_EQ(evaluate("sinh(0.5)"),  std::sinh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("cosh(0.5)"),  std::cosh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("tanh(0.5)"),  std::tanh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("asinh(0.5)"),  std::asinh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("acosh(2)"),    std::acosh(2));
  EXPECT_DOUBLE_EQ(evaluate("atanh(0.5)"),  std::atanh(0.5));
  EXPECT_DOUBLE_EQ(evaluate("erf(-1)"),   std::erf(-1));
  EXPECT_DOUBLE_EQ(evaluate("erfc(-1.5)"), std::erfc(-1.5));
  EXPECT_DOUBLE_EQ(evaluate("poltorectx(5, PI)"),    -5);
  EXPECT_DOUBLE_EQ(evaluate("poltorecty(5, PI)"),      5*std::sin(stk::expreval::pi()));
  EXPECT_DOUBLE_EQ(evaluate("recttopolr(-1, 0)"), 1);
  EXPECT_DOUBLE_EQ(evaluate("recttopola(-1, 0)"),  stk::expreval::pi());
  EXPECT_DOUBLE_EQ(evaluate("unit_step(0.5, 0, 1)"),   1);
  EXPECT_DOUBLE_EQ(evaluate("cycloidal_ramp(1, 0, 1)"),     1);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cos_ramp(1/3)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("cosine_ramp(1/3)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("linear_ramp(1/4, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("haversine_pulse(1/6, 0, 1)"),  0.25);
  EXPECT_DOUBLE_EQ(evaluate("point2d(1, 0, 1, 1)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("point3d(0, -1, 0, 1, 1)"),   0.5);
}

double reference_normal_pdf(double x, double mu, double sigma) {
  return std::exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma)) / std::sqrt(2.0*stk::expreval::pi()*sigma*sigma);
}

double reference_weibull_pdf(double x, double k, double lambda) {
  return (x >= 0) ? (k/lambda)*std::pow(x/lambda, k-1)*std::exp(-std::pow(x/lambda, k)) : 0;
}

double reference_gamma_pdf(double x, double k, double theta) {
  return (x >= 0) ? 1/(std::tgamma(k)*std::pow(theta, k))*std::pow(x, k-1)*std::exp(-x/theta) : 0;
}

TEST(HostEvaluation, testPDFFunctions)
{
  EXPECT_DOUBLE_EQ(evaluate("exponential_pdf(0, 1)"),  1);
  EXPECT_DOUBLE_EQ(evaluate("log_uniform_pdf(2, 1, E)"),   0.5);
  EXPECT_DOUBLE_EQ(evaluate("normal_pdf(0.75, 1, 0.5)"), reference_normal_pdf(0.75, 1, 0.5));
  EXPECT_DOUBLE_EQ(evaluate("weibull_pdf(1, 5, 1)"),    reference_weibull_pdf(1, 5, 1));
  EXPECT_DOUBLE_EQ(evaluate("gamma_pdf(5, 5, 1)"),   reference_gamma_pdf(5, 5, 1));
}
//-END

} // namespace <unnamed>

