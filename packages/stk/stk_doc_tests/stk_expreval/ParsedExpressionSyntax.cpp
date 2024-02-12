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
bool
isValidParse(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  EXPECT_NO_THROW(expr_eval.parse());
  return expr_eval.getSyntaxStatus();
}

bool
isInvalidParse(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  try {
    expr_eval.parse();
  }
  catch (std::runtime_error& ) {
    return !expr_eval.getSyntaxStatus();
  }

  return false;
}

TEST(ParsedEval, testAlgebraicSyntax)
{
  EXPECT_TRUE(isValidParse(""));
  EXPECT_TRUE(isValidParse(";;"));
  EXPECT_TRUE(isValidParse("2*2"));
  EXPECT_TRUE(isValidParse("3^2"));
  EXPECT_TRUE(isValidParse("x*-0.1"));
  EXPECT_TRUE(isValidParse("x*+0.1"));
  EXPECT_TRUE(isValidParse("x--7.0"));
  EXPECT_TRUE(isValidParse("x*-x"));
  EXPECT_TRUE(isValidParse("x*+x"));
  EXPECT_TRUE(isValidParse("v[0]=v[1]*0.1"));
  EXPECT_TRUE(isValidParse("x---x"));

  EXPECT_TRUE(isInvalidParse("0.01.02"));
  EXPECT_TRUE(isInvalidParse("5*.e+10"));
  EXPECT_TRUE(isInvalidParse("x y"));
  EXPECT_TRUE(isInvalidParse("x(y"));
  EXPECT_TRUE(isInvalidParse("x*"));
  EXPECT_TRUE(isInvalidParse("x*(y+1"));
  EXPECT_TRUE(isInvalidParse("cos(x"));
  EXPECT_TRUE(isInvalidParse("(x)y"));
  EXPECT_TRUE(isInvalidParse("()"));
}


bool
isValidFunction(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  EXPECT_NO_THROW(expr_eval.parse());
  return !expr_eval.undefinedFunction();
}

bool
isInvalidFunction(const char *expr)
{
  stk::expreval::Eval expr_eval(expr);
  try {
    expr_eval.parse();
  }
  catch (std::runtime_error& ) {
    return expr_eval.undefinedFunction();
  }

  return false;
}

TEST(ParsedEval, testFunctionSyntax)
{
  EXPECT_TRUE(isValidFunction("sin(1)"));
  EXPECT_TRUE(isValidFunction("SIN(1)"));
  EXPECT_TRUE(isValidFunction("rand()"));
  EXPECT_TRUE(isValidFunction("time()"));
  EXPECT_TRUE(isValidFunction("random()"));
  EXPECT_TRUE(isValidFunction("random(1)"));
  EXPECT_TRUE(isValidFunction("cosine_ramp(x,y)"));
  EXPECT_TRUE(isValidFunction("normal_pdf(x, alpha, beta)"));

  EXPECT_TRUE(isInvalidFunction("stress(1)"));
  EXPECT_TRUE(isInvalidFunction("gamma(1)"));
}
//-END
} // namespace <unnamed>
