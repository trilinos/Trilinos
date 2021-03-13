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
#include <stk_util/parallel/ParallelComm.hpp>

namespace {

class TestCommBuffer : public ::testing::Test
{
public:
  template<typename PACK_LAMBDA>
  void pack(const PACK_LAMBDA& pack_stuff)
  {
    pack_stuff(buf);
    allocate_buf();
    pack_stuff(buf);
  }

  void allocate_buf()
  {
    data.resize(buf.size());
    buf.set_buffer_ptrs(data.data(), data.data(), data.data()+data.size());
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack(const T& gold)
  {
    T val;
    buf.unpack(val);
    EXPECT_NEAR(gold, val, 1.e-6);
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value>::type
  check_unpack(const std::vector<T>& gold)
  {
    std::vector<T> val;
    buf.unpack(val);
    ASSERT_EQ(val.size(), gold.size());
    for(unsigned i=0; i<gold.size(); ++i) {
      EXPECT_NEAR(gold[i], val[i], 1.e-6);
    }
  }

  template<typename T>
  typename std::enable_if<!std::is_floating_point<T>::value>::type
  check_unpack(const T& gold)
  {
    T val;
    buf.unpack(val);
    EXPECT_EQ(gold, val);
  }

  template<typename T>
  void check_peek_and_unpack_pair(const T& gold)
  {
    T peekVal;
    buf.peek<T>(peekVal);
    EXPECT_EQ(gold, peekVal);
    T val;
    buf.unpack<T>(val);
    EXPECT_EQ(gold, val);
  }

  stk::CommBuffer buf;
  std::vector<unsigned char> data;
};

TEST_F(TestCommBuffer, pack_unpack_primitives)
{
  char goldc = 's';
  int goldn = 5;
  float goldf = 1.1;
  double goldd = 3.3;
  bool goldb = true;

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldc);
    buffer.pack(goldn);
    buffer.pack(goldf);
    buffer.pack(goldd);
    buffer.pack(goldb);
  });

  buf.reset();

  check_unpack(goldc);
  check_unpack(goldn);
  check_unpack(goldf);
  check_unpack(goldd);
  check_unpack(goldb);
}

TEST_F(TestCommBuffer, pack_unpack_pairs)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpic);
    buffer.pack(goldpdi);
  });

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, skip_pack_pairs_explicit_template_parameter)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  buf.skip<std::pair<int,char>>(1);
  buf.skip<std::pair<unsigned long,int>>(1);
  buf.skip<std::pair<double,int>>(1);
  
  allocate_buf();

  buf.pack<std::pair<int,char>>(goldpic);
  buf.pack<std::pair<unsigned long,int>>(goldpui);
  buf.pack<std::pair<double,int>>(goldpdi);

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpui);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pairs_explicit_template_parameter)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack<std::pair<int,char>>(goldpic);
    buffer.pack<std::pair<unsigned long,int>>(goldpui);
    buffer.pack<std::pair<double,int>>(goldpdi);
  });

  buf.reset();

  check_unpack(goldpic);
  check_unpack(goldpui);
  check_unpack(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pairs_explicit_template_parameter_peek_unpack)
{
  std::pair<int,char> goldpic(5, 'c');
  std::pair<unsigned long,int> goldpui(5, 9);
  std::pair<double,int> goldpdi(3.14, 42);

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpic);
    buffer.pack(goldpui);
    buffer.pack(goldpdi);
  });

  buf.reset();

  check_peek_and_unpack_pair(goldpic);
  check_peek_and_unpack_pair(goldpui);
  check_peek_and_unpack_pair(goldpdi);
}

TEST_F(TestCommBuffer, pack_unpack_pair_string)
{
  std::pair<std::string,int> goldpsi("short string", 5);
  std::pair<double,std::string> goldpds(3.14, "some string");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldpsi);
    buffer.pack(goldpds);
  });

  buf.reset();

  check_unpack(goldpsi);
  check_unpack(goldpds);
}

TEST_F(TestCommBuffer, pack_unpack_vector_primitives)
{
  std::vector<char> goldc = {'a', 'b', 'c'};
  std::vector<int> goldn = {5, 6, 7};
  std::vector<float> goldf = {1.1, 1.2, 1.3};
  std::vector<double> goldd = {3.3, 3.4, 3.5};
  std::vector<bool> goldb = {true, false, true};

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldc);
    buffer.pack(goldn);
    buffer.pack(goldf);
    buffer.pack(goldd);
    buffer.pack(goldb);
  });

  buf.reset();

  check_unpack(goldc);
  check_unpack(goldn);
  check_unpack(goldf);
  check_unpack(goldd);
  check_unpack(goldb);
}

TEST_F(TestCommBuffer, pack_unpack_string)
{
  std::string goldstr("stk is awesome");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldstr);
  });

  buf.reset();

  check_unpack(goldstr);
}

TEST_F(TestCommBuffer, pack_unpack_string_explicit_template_parameter)
{
  std::string goldstr("stk is awesome");

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack<std::string>(goldstr);
  });

  buf.reset();

  std::string str;
  buf.unpack<std::string>(str);
  EXPECT_EQ(goldstr, str);
}

TEST_F(TestCommBuffer, pack_unpack_vector_string)
{
  std::vector<std::string> goldstr = {"stk", "is", "awesome"};

  pack([=](stk::CommBuffer& buffer) {
    buffer.pack(goldstr);
  });

  buf.reset();

  check_unpack(goldstr);
}

}

