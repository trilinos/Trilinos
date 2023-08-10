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
TEST(HostEvaluation, bindScalar)
{
  stk::expreval::Eval expr("x=5; y=y+x; y+z");
  expr.parse();
  double y = 3.0;
  double z = 4.0;
  expr.bindVariable("y", y, 1);
  expr.bindVariable("z", z, 1);
  EXPECT_DOUBLE_EQ(expr.evaluate(), 12);
}

TEST(HostEvaluation, bindVector)
{
  stk::expreval::Eval expr("(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])^0.5");
  expr.parse();
  double a[3] = {1, 2, 3};
  double b[3] = {5, 4, 4};
  expr.bindVariable("a", *a, 3);
  expr.bindVariable("b", *b, 3);
  EXPECT_DOUBLE_EQ(expr.evaluate(), 5);
}

TEST(HostEvaluation, bindVectorOneBasedIndex)
{
  stk::expreval::Eval expr("(a[1]*b[1] + a[2]*b[2] + a[3]*b[3])^0.5", stk::expreval::Variable::ONE_BASED_INDEX);
  expr.parse();
  double a[3] = {1, 2, 3};
  double b[3] = {5, 4, 4};
  expr.bindVariable("a", *a, 3);
  expr.bindVariable("b", *b, 3);
  EXPECT_DOUBLE_EQ(expr.evaluate(), 5);
}
//-END

} // namespace <unnamed>
