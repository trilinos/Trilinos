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

#include "gtest/gtest.h"                // for ASSERT_TRUE, TEST
#include "stk_mesh/baseImpl/DoOp.hpp"
#include <limits>
#include <iomanip>

TEST(UnitTestDoOp, basic_sum)
{
  EXPECT_FLOAT_EQ(5.9, (stk::mesh::impl::DoOp<float,stk::mesh::Operation::SUM>()(2.2, 3.7)));
  EXPECT_DOUBLE_EQ(5.9, (stk::mesh::impl::DoOp<double,stk::mesh::Operation::SUM>()(2.2, 3.7)));
  EXPECT_EQ(5, (stk::mesh::impl::DoOp<int,stk::mesh::Operation::SUM>()(2, 3)));
}

TEST(UnitTestDoOp, basic_min)
{
  EXPECT_FLOAT_EQ(2.2, (stk::mesh::impl::DoOp<float,stk::mesh::Operation::MIN>()(2.2, 3.7)));
  EXPECT_DOUBLE_EQ(2.2, (stk::mesh::impl::DoOp<double,stk::mesh::Operation::MIN>()(2.2, 3.7)));
  EXPECT_EQ(2, (stk::mesh::impl::DoOp<int,stk::mesh::Operation::MIN>()(2, 3)));
}

TEST(UnitTestDoOp, basic_max)
{
  EXPECT_FLOAT_EQ(3.7, (stk::mesh::impl::DoOp<float,stk::mesh::Operation::MAX>()(2.2, 3.7)));
  EXPECT_DOUBLE_EQ(3.7, (stk::mesh::impl::DoOp<double,stk::mesh::Operation::MAX>()(2.2, 3.7)));
  EXPECT_EQ(3, (stk::mesh::impl::DoOp<int,stk::mesh::Operation::MAX>()(2, 3)));
}

TEST(UnitTestDoOp, sum_large_small)
{
  double small = -2.0518531544667047e-02;
  double large = 2.5372303414487405e+03;
  
  double sum = stk::mesh::impl::DoOp<double,stk::mesh::Operation::SUM>()(large, small);

  const double expectedSum = 2.5372098229171961e+03;
  const double error = std::abs(expectedSum - sum);
  constexpr double expectedError = 2.e-13;

  EXPECT_TRUE(error <= expectedError);
}
