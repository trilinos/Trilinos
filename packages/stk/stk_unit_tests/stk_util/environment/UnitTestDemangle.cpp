// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stk_util/environment/Demangle.hpp>  // for demangle
#include <gtest/gtest.h>
#include <string>                       // for operator==, string, etc
#include <typeinfo>                     // for type_info
#include <vector>                       // for vector



bool
utest_demangle()
{
  return true;
}

TEST(UnitTestDemangle, UnitTest)
{
#if defined(__PGI)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool ()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    ASSERT_EQ(linux_name, demangled_name);
  }

#elif defined(__sun)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    ASSERT_EQ(linux_name, demangled_name);
  }

#elif defined(__xlC__)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
//    ASSERT_EQ((linux_name == demangled_name), true);
  }

  {
    std::string linux_name("bool ()()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
//    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double> >");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
//    ASSERT_EQ((linux_name == demangled_name), true);
  }
#elif defined(__linux__)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    ASSERT_EQ((linux_name == demangled_name), true);
  }

  {
#if (__GNUC_MINOR__ > 5) || defined(__clang__)
    std::string linux_name("bool ()");
#else
    std::string linux_name("bool ()()");
#endif
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    ASSERT_EQ(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

#ifdef _GLIBCXX_DEBUG
    std::string linux_name("__gnu_debug_def::vector<double, std::allocator<double> >");
#else
    std::string linux_name("std::vector<double, std::allocator<double> >");
#endif
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    ASSERT_EQ((linux_name == demangled_name), true);
  }
#endif
}
