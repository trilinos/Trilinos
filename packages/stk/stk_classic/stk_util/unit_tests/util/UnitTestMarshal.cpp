/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>

#include <stk_util/util/Marshal.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

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
stk_classic::Marshal &operator<<(stk_classic::Marshal &mout, const S1 &s) {
  return mout << s.m_string << s.m_int;
}

stk_classic::Marshal &operator>>(stk_classic::Marshal &min, S1 &s) {
  return min >> s.m_string >> s.m_int;
}


// Define Marshal for S2 type in stk's namespace
namespace stk_classic {

template<>
stk_classic::Marshal &operator<<(stk_classic::Marshal &mout, const S2 &s) {
  return mout << s.m_string << s.m_int;
}

template<>
stk_classic::Marshal &operator>>(stk_classic::Marshal &min, S2 &s) {
  return min >> s.m_string >> s.m_int;
}

} // namespace stk_classic


template <class T>
void
test(const T &t_in) 
{
  stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_ALL);
  mout << t_in;

  stk_classic::Marshal min(mout.str());
  T t_out;
  min >> t_out;

  STKUNIT_ASSERT_EQUAL(t_in, t_out);
}

STKUNIT_UNIT_TEST(UnitTestMarshal, UnitTest)
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
    stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_ALL);
    std::string s_out("this is a test");
    int i_out = 7;
      
    mout << s_out << i_out;

    stk_classic::Marshal min(mout.str());
    std::string s_in;
    int i_in;

    min >> s_in >> i_in;

    STKUNIT_ASSERT_EQUAL(min.size(), mout.size());
    STKUNIT_ASSERT(mout);
    STKUNIT_ASSERT(min);
    STKUNIT_ASSERT_EQUAL((s_in == s_out), true);
    STKUNIT_ASSERT_EQUAL(i_in, i_out);
  }

  // Marshal/Unmarshal locally defined class/struct
  {
    stk_classic::Marshal mout;
    S1 s_out("this is a test", 5);
    mout << s_out;

    stk_classic::Marshal min(mout.str());
    S1 s_in;

    min >> s_in;

    STKUNIT_ASSERT_EQUAL((s_in.m_string == s_out.m_string), true);
    STKUNIT_ASSERT_EQUAL(s_in.m_int, s_out.m_int);
  }

  // Marshal/Unmarshal std::vector of S2's (uses the stk namespace operator>> and operator<<)
  {
    stk_classic::Marshal mout;
    S2 s_out("this is a test", 5);
    std::vector<S2> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk_classic::Marshal min(mout.str());
    std::vector<S2> v_in;
    
    min >> v_in;

    STKUNIT_ASSERT_EQUAL((v_in[0].m_string == v_out[0].m_string), true);
  }

  // Marshal/Unmarshal error from type mismatch
  {
    
    stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_ALL);
    std::string s_out("this is a test");
    int i_out = 7;
      
    mout << s_out << i_out;

    stk_classic::Marshal min(mout.str());
    std::string s_in;
    double x_in;

    STKUNIT_ASSERT_THROW(min >> s_in >> x_in, std::runtime_error);
  }
  
  // Marshal error for STL container mismatch
  {
    stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_ALL);
    S1 s_out("this is a test", 5);
    std::vector<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk_classic::Marshal min(mout.str());
    std::list<S1> v_in;
    
    STKUNIT_ASSERT_THROW(min >> v_in, std::runtime_error);
  }
  
  // Marshal error for STL container mismatch
  {
    stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_ALL);
    S1 s_out("this is a test", 5);
    std::list<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk_classic::Marshal min(mout.str());
    std::vector<S1> v_in;
    
    STKUNIT_ASSERT_THROW(min >> v_in, std::runtime_error);
  }
  
  // Marshal without error for STL container mismatch
  {
    stk_classic::Marshal mout(stk_classic::Marshal::TYPE_CHECK_NONE);
    S1 s_out("this is a test", 5);
    std::vector<S1> v_out;
    v_out.push_back(s_out);
    
    mout << v_out;

    stk_classic::Marshal min(mout.str());
    std::list<S1> v_in;
    
    STKUNIT_ASSERT_NO_THROW(min >> v_in);
  }
}
