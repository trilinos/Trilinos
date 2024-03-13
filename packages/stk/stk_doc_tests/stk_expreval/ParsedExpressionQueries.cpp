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
TEST(ParsedEval, isConstantExpression)
{
  stk::expreval::Eval evalEmpty;
  evalEmpty.parse();
  EXPECT_TRUE(evalEmpty.is_constant_expression());

  stk::expreval::Eval evalConstant("2");
  evalConstant.parse();
  EXPECT_TRUE(evalConstant.is_constant_expression());

  stk::expreval::Eval evalVar("x");
  evalVar.parse();
  EXPECT_FALSE(evalVar.is_constant_expression());
}

TEST(ParsedEval, isVariable)
{
  stk::expreval::Eval evalEmpty;
  evalEmpty.parse();
  EXPECT_FALSE(evalEmpty.is_variable("x"));
  
  stk::expreval::Eval evalTwoVar("x + y");
  evalTwoVar.parse();
  EXPECT_TRUE(evalTwoVar.is_variable("x"));
  EXPECT_TRUE(evalTwoVar.is_variable("y"));
  EXPECT_FALSE(evalTwoVar.is_variable("z"));

  stk::expreval::Eval evalInsVar("lambda + Lambda");
  evalInsVar.parse();
  EXPECT_EQ(evalInsVar.get_variable_names().size(), 1u);
  EXPECT_TRUE(evalInsVar.is_variable("LAMBDA"));
  EXPECT_TRUE(evalInsVar.is_variable("lambda"));
  EXPECT_TRUE(evalInsVar.is_variable("Lambda"));

}

TEST(ParsedEval, isScalar)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_TRUE(eval.is_scalar("x"));

  stk::expreval::Eval evalBind("y^2");
  evalBind.parse();
  EXPECT_TRUE(evalBind.is_scalar("y"));

  stk::expreval::Eval evalBindArray("z");
  evalBindArray.parse();
  double z[3] = {4.0, 5.0, 6.0};
  evalBindArray.bindVariable("z", *z, 3);
  EXPECT_FALSE(evalBindArray.is_scalar("z"));
}

TEST(ParsedEval, getAllVariables)
{
  stk::expreval::Eval eval;
  eval.parse();
  EXPECT_EQ(eval.get_variable_names().size(), 0u);

  stk::expreval::Eval evalVars("x = sin(y)");
  evalVars.parse();
  EXPECT_EQ(evalVars.get_variable_names().size(), 2u);
  EXPECT_TRUE(evalVars.is_variable("x"));
  EXPECT_TRUE(evalVars.is_variable("y"));
}

TEST(ParsedEval, getDependentVariables)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_EQ(eval.get_dependent_variable_names().size(), 0u);

  stk::expreval::Eval evalAssign("x = 2");
  evalAssign.parse();
  EXPECT_EQ(evalAssign.get_dependent_variable_names().size(), 1u);
  EXPECT_TRUE(evalAssign.is_variable("x"));

  stk::expreval::Eval evalTwoVar("x = 2; y = x");
  evalTwoVar.parse();
  EXPECT_EQ(evalTwoVar.get_dependent_variable_names().size(), 2u);
  EXPECT_TRUE(evalTwoVar.is_variable("x"));
  EXPECT_TRUE(evalTwoVar.is_variable("y"));
}

TEST(ParsedEval, getIndependentVariables)
{
  stk::expreval::Eval eval("x");
  eval.parse();
  EXPECT_EQ(eval.get_independent_variable_names().size(), 1u);
  EXPECT_TRUE(eval.is_variable("x"));

  stk::expreval::Eval evalAssign("x = 2");
  evalAssign.parse();
  EXPECT_EQ(evalAssign.get_independent_variable_names().size(), 0u);
  EXPECT_TRUE(evalAssign.is_variable("x"));

  stk::expreval::Eval evalTwoVar("x = sin(y)");
  evalTwoVar.parse();
  EXPECT_EQ(evalTwoVar.get_independent_variable_names().size(), 1u);
  EXPECT_TRUE(evalTwoVar.is_variable("x"));
  EXPECT_TRUE(evalTwoVar.is_variable("y"));
}

//-END

} // namespace <unnamed>
