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
TEST(Variable, scalar_vs_array)
{
  stk::expreval::Eval expr("x[1]", stk::expreval::Variable::ArrayOffset::ZERO_BASED_INDEX);
  expr.parse();
  EXPECT_TRUE(expr.is_scalar("x"));
  EXPECT_EQ(expr.getValue("x"), 0.0);
  EXPECT_ANY_THROW(expr.evaluate());

  double x[2] = {3.0, 4.0};
  expr.bindVariable("x", *x, 2);
  EXPECT_FALSE(expr.is_scalar("x"));
  EXPECT_EQ(expr.evaluate(), 4.0);

  stk::expreval::Eval expr2("y[1]", stk::expreval::Variable::ArrayOffset::ONE_BASED_INDEX);
  expr2.parse();
  double y[2] = {3.0, 4.0};
  expr2.bindVariable("y", *y, 2);
  EXPECT_FALSE(expr2.is_scalar("y"));
  EXPECT_EQ(expr2.evaluate(), 3.0);
}

TEST(Variable, demonstrate_states)
{
  stk::expreval::Eval expr("x");
  expr.parse();
  EXPECT_EQ(expr.evaluate(), 0.0);

  double x = 2.0;
  expr.bindVariable("x", x, 1);
  EXPECT_EQ(expr.evaluate(), 2.0);

  expr.unbindVariable("x");
  EXPECT_EQ(expr.evaluate(), 0.0);

  expr.deactivateVariable("x");
  EXPECT_ANY_THROW(expr.evaluate());
}
//-END
} // namespace <unnamed>
