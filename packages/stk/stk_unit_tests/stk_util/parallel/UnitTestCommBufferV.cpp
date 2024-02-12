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

#include "gtest/gtest.h"
#include "stk_util/environment/WallTime.hpp"  // for wall_time
#include "stk_util/parallel/CommBufferV.hpp"  // for CommBufferV
#include <cstddef>                            // for size_t
#include <iostream>                           // for operator<<, endl, basic_ostream, basic_ostr...
#include <vector>                             // for vector


TEST(Parallel, comm_buffer_exp_double)
{
  const size_t num_doubles = 8;
  double data1 = 99.0;
  std::vector<double> doubles(num_doubles, data1);
  const double expected_data1 = 99.0;

  stk::CommBufferV buf;

  buf.pack<double>(data1);
  size_t expected_size_in_bytes = sizeof(double);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());

  buf.pack<double>(doubles.data(), num_doubles);
  expected_size_in_bytes += num_doubles*sizeof(double);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());

  buf.unpack<double>(data1);
  expected_size_in_bytes -= sizeof(double);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());
  EXPECT_EQ(expected_data1, data1);

  buf.unpack<double>(doubles.data(), num_doubles);
  expected_size_in_bytes = 0;
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());
  for(double data : doubles) {
    EXPECT_EQ(expected_data1, data);
  }
}

TEST(Parallel, comm_buffer_exp_mixed)
{
  double double1 = 99.9;
  int int1 = 9;
  const size_t num = 8;
  std::vector<double> doubles(num, double1);
  std::vector<int> ints(num, int1);
  const double expected_double1 = 99.9;
  const int expected_int1 = 9;

  stk::CommBufferV buf;

  buf.pack<double>(double1);
  buf.pack<int>(int1);
  size_t expected_size_in_bytes = sizeof(double)+sizeof(int);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());

  buf.pack<double>(doubles.data(), num);
  buf.pack<int>(ints.data(), num);
  expected_size_in_bytes += num*(sizeof(double)+sizeof(int));
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());

  buf.unpack<double>(double1);
  expected_size_in_bytes -= sizeof(double);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());
  EXPECT_EQ(expected_double1, double1);

  buf.unpack<int>(int1);
  expected_size_in_bytes -= sizeof(int);
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());
  EXPECT_EQ(expected_int1, int1);

  buf.unpack<double>(doubles.data(), num);
  buf.unpack<int>(ints.data(), num);
  expected_size_in_bytes = 0;
  EXPECT_EQ(expected_size_in_bytes, buf.size_in_bytes());
  for(double data : doubles) {
    EXPECT_EQ(expected_double1, data);
  }
  for(double intdata : ints) {
    EXPECT_EQ(expected_int1, intdata);
  }
}

