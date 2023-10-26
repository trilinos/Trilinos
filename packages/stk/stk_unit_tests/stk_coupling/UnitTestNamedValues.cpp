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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_coupling/impl_NamedValues.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <cstdint>
#include <string>
#include <vector>
#include <map>

namespace {

TEST(UnitTestNamedValues, get_and_set_and_has)
{
  stk::coupling::impl::NamedValues values;

  bool boolValue = true;
  int intValue = 3;
  float floatValue = 1.14;
  double doubleValue = 3.14;
  std::string stringValue("foo");
  std::string stringLiteralValue("hmmm");
  std::vector<int> vectorIntValues = {1, 2, 3};
  std::vector<float> vectorFloatValues = {0.1, 0.2, 0.3};
  std::vector<double> vectorDoubleValues = {1.1, 2.2, 3.3};
  std::vector<std::string> vectorStringValues = {"one", "two", "three"};
  std::vector<std::pair<std::string,int>> vectorPairStringIntValues =
   {{"one",1}, {"two",2}, {"three",3}};

  values.set_value("boolValue", boolValue);
  values.set_value("intValue", intValue);
  values.set_value("floatValue", floatValue);
  values.set_value("doubleValue", doubleValue);
  values.set_value("stringValue", stringValue);
  //values.set_value("stringLiteralValue", "hmmm"); //won't compile
  values.set_value<std::string>("stringLiteralValue", "hmmm");
  values.set_value("vectorIntValue", vectorIntValues);
  values.set_value("vectorFloatValue", vectorFloatValues);
  values.set_value("vectorDoubleValue", vectorDoubleValues);
  values.set_value("vectorStringValue", vectorStringValues);
  values.set_value("vectorPairStringIntValue", vectorPairStringIntValues);

  EXPECT_FALSE(values.has_value<bool>("garbage"));
  EXPECT_FALSE(values.has_value<int>("garbage"));
  EXPECT_FALSE(values.has_value<float>("garbage"));
  EXPECT_FALSE(values.has_value<double>("garbage"));
  EXPECT_FALSE(values.has_value<std::string>("garbage"));
  EXPECT_FALSE(values.has_value<std::vector<int>>("garbage"));
  EXPECT_FALSE(values.has_value<std::vector<float>>("garbage"));
  EXPECT_FALSE(values.has_value<std::vector<double>>("garbage"));
  EXPECT_FALSE(values.has_value<std::vector<std::string>>("garbage"));
  EXPECT_FALSE((values.has_value<std::vector<std::pair<std::string,int>>>("garbage")));

  EXPECT_TRUE(values.has_value<bool>("boolValue"));
  EXPECT_TRUE(values.has_value<int>("intValue"));
  EXPECT_TRUE(values.has_value<float>("floatValue"));
  EXPECT_TRUE(values.has_value<double>("doubleValue"));
  EXPECT_TRUE(values.has_value<std::string>("stringValue"));
  EXPECT_FALSE(values.has_value<char*>("stringLiteralValue"));
  EXPECT_TRUE(values.has_value<std::string>("stringLiteralValue"));
  EXPECT_TRUE(values.has_value<std::vector<int>>("vectorIntValue"));
  EXPECT_TRUE(values.has_value<std::vector<float>>("vectorFloatValue"));
  EXPECT_TRUE(values.has_value<std::vector<double>>("vectorDoubleValue"));
  EXPECT_TRUE(values.has_value<std::vector<std::string>>("vectorStringValue"));
  EXPECT_TRUE((values.has_value<std::vector<std::pair<std::string,int>>>("vectorPairStringIntValue")));

  EXPECT_ANY_THROW(values.get_value<bool>("garbage"));
  EXPECT_ANY_THROW(values.get_value<int>("garbage"));
  EXPECT_ANY_THROW(values.get_value<float>("garbage"));
  EXPECT_ANY_THROW(values.get_value<double>("garbage"));
  EXPECT_ANY_THROW(values.get_value<std::string>("garbage"));
  EXPECT_ANY_THROW(values.get_value<std::vector<int>>("garbage"));
  EXPECT_ANY_THROW(values.get_value<std::vector<float>>("garbage"));
  EXPECT_ANY_THROW(values.get_value<std::vector<double>>("garbage"));
  EXPECT_ANY_THROW(values.get_value<std::vector<std::string>>("garbage"));
  EXPECT_ANY_THROW((values.get_value<std::vector<std::pair<std::string,int>>>("garbage")));

  EXPECT_EQ(boolValue, values.get_value<bool>("boolValue"));
  EXPECT_EQ(intValue, values.get_value<int>("intValue"));
  EXPECT_NEAR(floatValue, values.get_value<float>("floatValue"), 1.e-6);
  EXPECT_NEAR(doubleValue, values.get_value<double>("doubleValue"), 1.e-12);
  EXPECT_EQ(stringValue, values.get_value<std::string>("stringValue"));
  EXPECT_EQ(stringLiteralValue, values.get_value<std::string>("stringLiteralValue"));
  EXPECT_EQ(vectorIntValues, values.get_value<std::vector<int>>("vectorIntValue"));
  EXPECT_EQ(vectorFloatValues, values.get_value<std::vector<float>>("vectorFloatValue"));
  EXPECT_EQ(vectorDoubleValues, values.get_value<std::vector<double>>("vectorDoubleValue"));
  EXPECT_EQ(vectorStringValues, values.get_value<std::vector<std::string>>("vectorStringValue"));
  EXPECT_EQ(vectorPairStringIntValues, (values.get_value<std::vector<std::pair<std::string,int>>>("vectorPairStringIntValue")));
}

TEST(UnitTestNamedValues, comparison)
{
  stk::coupling::impl::NamedValues values;

  bool boolValue = true;
  int intValue = 3;
  double floatValue = 1.14;
  double doubleValue = 3.14;
  std::string stringValue("foo");
  std::vector<float> vectorFloatValue = {0.3, 0.4, 0.5};

  values.set_value("boolValue", boolValue);
  values.set_value("intValue", intValue);
  values.set_value("floatValue", floatValue);
  values.set_value("doubleValue", doubleValue);
  values.set_value("stringValue", stringValue);
  values.set_value("vectorFloatValue", vectorFloatValue);

  stk::coupling::impl::NamedValues otherValues;

  EXPECT_FALSE(otherValues == values);

  otherValues.set_value("boolValue", boolValue);
  otherValues.set_value("intValue", 4);
  otherValues.set_value("floatValue", floatValue);
  otherValues.set_value("doubleValue", doubleValue);
  otherValues.set_value("stringValue", stringValue);
  otherValues.set_value("vectorFloatValue", vectorFloatValue);

  EXPECT_FALSE(otherValues == values);

  otherValues.set_value("intValue", intValue);

  EXPECT_TRUE(otherValues == values);
}

struct TestCommBuffer
{

  void allocate()
  {
    data.resize(buffer.size());
    buffer.set_buffer_ptrs(data.data(), data.data(), data.data()+data.size());
  }

  stk::CommBuffer buffer;
  std::vector<unsigned char> data;
};

TEST(UnitTestNamedValues, pack_unpack)
{
  stk::coupling::impl::NamedValues values;

  bool boolValue = true;
  int intValue = 3;
  float floatValue = 1.14;
  double doubleValue = 3.14;
  std::string stringValue("foo");
  std::vector<int> vectorIntValues = {3, 2, 1};
  std::vector<float> vectorFloatValues = {1.1, 2.2, 3.3};
  std::vector<double> vectorDoubleValues = {1.1, 2.2, 3.3};
  std::vector<std::string> vectorStringValues = {"one", "two", "three"};
  std::vector<std::pair<std::string,int>> vectorPairStringIntValues =
   {{"one",1}, {"two",2}, {"three",3}};

  values.set_value("boolValue", boolValue);
  values.set_value("intValue", intValue);
  values.set_value("floatValue", floatValue);
  values.set_value("doubleValue", doubleValue);
  values.set_value("stringValue", stringValue);
  values.set_value("vectorIntValue", vectorIntValues);
  values.set_value("vectorFloatValue", vectorFloatValues);
  values.set_value("vectorDoubleValue", vectorDoubleValues);
  values.set_value("vectorStringValue", vectorStringValues);
  values.set_value("vectorPairStringIntValue", vectorPairStringIntValues);

  TestCommBuffer testCommBuf;
  stk::CommBuffer& buf = testCommBuf.buffer;

  values.pack(buf);
  testCommBuf.allocate();
  values.pack(buf);
  buf.reset(); //resets internal pointers to prepare for unpacking.

  stk::coupling::impl::NamedValues otherValues;
  otherValues.unpack(buf);

  EXPECT_TRUE(values == otherValues);
}

TEST(UnitTestNamedValues, commSparse)
{
  std::vector<std::pair<std::string,int>> vectorPairStringIntValues =
   {{"one",1}, {"two",2}, {"three",3}};

  stk::coupling::impl::NamedValues values;
  values.set_value("vectorPairStringIntValue", vectorPairStringIntValues);

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::CommSparse commSparse(comm);
  int numProcs = stk::parallel_machine_size(comm);
  int myRank = stk::parallel_machine_rank(comm);
  int destRank = (myRank+1) % numProcs;

  const bool needToUnpackRecvdMessage =
    stk::pack_and_communicate(commSparse, [&commSparse, &destRank, &values]() {
      values.pack(commSparse.send_buffer(destRank));
    });

  EXPECT_TRUE(needToUnpackRecvdMessage);
  stk::coupling::impl::NamedValues otherValues;
  int sourceRank = (myRank-1+numProcs) % numProcs;
  otherValues.unpack(commSparse.recv_buffer(sourceRank));
  EXPECT_EQ(values, otherValues);
}

}
