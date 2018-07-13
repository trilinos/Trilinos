// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <list>                         // for list
#include <stdexcept>                    // for runtime_error
#include <gtest/gtest.h>
#include <stk_util/util/Marshal.hpp>    // for Marshal, operator<<, etc
#include <string>                       // for operator==, string, etc
#include <vector>                       // for vector



struct S1
{
  S1()
    : m_string(),
      m_int(0)
  {}
  
  S1(const std::string &s, int i)
    : m_string(s),
      m_int(i)
  {}
  
  std::string   m_string;
  int           m_int;
};

struct S2
{
  S2()
    : m_string(),
      m_int(0)
  {}
  
  S2(const std::string &s, int i)
    : m_string(s),
      m_int(i)
  {}
  
  std::string   m_string;
  int           m_int;
};

// Define Marshal for S1 type in S1's namespace
stk::Marshal &operator<<(stk::Marshal &mout, const S1 &s) {
  return mout << s.m_string << s.m_int;
}

stk::Marshal &operator>>(stk::Marshal &min, S1 &s) {
  return min >> s.m_string >> s.m_int;
}


// Define Marshal for S2 type in stk's namespace
namespace stk {

template<>
stk::Marshal &operator<<(stk::Marshal &mout, const S2 &s) {
  return mout << s.m_string << s.m_int;
}

template<>
stk::Marshal &operator>>(stk::Marshal &min, S2 &s) {
  return min >> s.m_string >> s.m_int;
}

} // namespace stk


template <class T>
void
test(const T &t_in) 
{
  stk::Marshal mout(stk::Marshal::TYPE_CHECK_ALL);
  mout << t_in;

  stk::Marshal min(mout.str());
  T t_out;
  min >> t_out;

  ASSERT_EQ(t_in, t_out);
}

TEST(UnitTestMarshal, UnitTest)
{
  // Marshal/Unmarshal POD
  {
    test<char>('a');
    test<signed char>('a');
    test<unsigned char>('a');
    test<signed short>(1);
    test<unsigned short>(2);
    test<signed int>(3);
    test<unsigned int>(4);
    test<signed long>(5);
    test<unsigned long>(6);
    test<signed long long>(7);
    test<unsigned long long>(8);
    test<float>(9.0);
    test<double>(10.0);
  }

  // Marshal/Unmarshal more than one POD in a single message
  {
    stk::Marshal mout(stk::Marshal::TYPE_CHECK_ALL);
    std::string s_out("this is a test");
    int i_out = 7;
      
    mout << s_out << i_out;

    stk::Marshal min(mout.str());
    std::string s_in;
    int i_in = 0;

    min >> s_in >> i_in;

    ASSERT_EQ(min.size(), mout.size());
    ASSERT_TRUE(mout);
    ASSERT_TRUE(min);
    ASSERT_EQ((s_in == s_out), true);
    ASSERT_EQ(i_in, i_out);
  }

  // Marshal/Unmarshal locally defined class/struct
  {
    stk::Marshal mout;
    S1 s_out("this is a test", 5);
    mout << s_out;

    stk::Marshal min(mout.str());
    S1 s_in;

    min >> s_in;

    ASSERT_EQ((s_in.m_string == s_out.m_string), true);
    ASSERT_EQ(s_in.m_int, s_out.m_int);
  }

  // Marshal/Unmarshal std::vector of S2's (uses the stk namespace operator>> and operator<<)
  {
    stk::Marshal mout;
    S2 s_out("this is a test", 5);
    std::vector<S2> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk::Marshal min(mout.str());
    std::vector<S2> v_in;
    
    min >> v_in;

    ASSERT_EQ((v_in[0].m_string == v_out[0].m_string), true);
  }

  // Marshal/Unmarshal error from type mismatch
  {
    
    stk::Marshal mout(stk::Marshal::TYPE_CHECK_ALL);
    std::string s_out("this is a test");
    int i_out = 7;
      
    mout << s_out << i_out;

    stk::Marshal min(mout.str());
    std::string s_in;
    double x_in = 0.0;

    ASSERT_THROW(min >> s_in >> x_in, std::runtime_error);
  }
  
  // Marshal error for STL container mismatch
  {
    stk::Marshal mout(stk::Marshal::TYPE_CHECK_ALL);
    S1 s_out("this is a test", 5);
    std::vector<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk::Marshal min(mout.str());
    std::list<S1> v_in;
    
    ASSERT_THROW(min >> v_in, std::runtime_error);
  }
  
  // Marshal error for STL container mismatch
  {
    stk::Marshal mout(stk::Marshal::TYPE_CHECK_ALL);
    S1 s_out("this is a test", 5);
    std::list<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk::Marshal min(mout.str());
    std::vector<S1> v_in;
    
    ASSERT_THROW(min >> v_in, std::runtime_error);
  }
  
  // Marshal without error for STL container mismatch
  {
    stk::Marshal mout(stk::Marshal::TYPE_CHECK_NONE);
    S1 s_out("this is a test", 5);
    std::vector<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk::Marshal min(mout.str());
    std::list<S1> v_in;
    
    ASSERT_NO_THROW(min >> v_in);
  }
}
