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

#include <stk_util/util/ci_string.hpp>


std::ostream &
operator<<(
  std::ostream &        os,
  const ci_string &      s)
{
  return os << std::string(s.data(),s.length());
}
 

std::istream &
operator>>(
  std::istream &        is,
  ci_string &            s)
{
  std::string t;
  is >> t;
  s = ci_string(t.begin(), t.end());
  
  return is;
}

std::string
operator+(
  const std::string &   s1,
  const ci_string &      s2) 
{
  std::string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


ci_string
operator+(
  const char *          s1,
  const ci_string &      s2) 
{
  ci_string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


ci_string
operator+(
  const ci_string &      s1,
  const std::string &   s2) 
{
  ci_string s(s1);
  
  return s.append(s2.begin(), s2.end());
}


// namespace boost {

// std::size_t
// hash<ci_string>::operator()(
//   const ci_string &      s) const
// {
//   std::size_t seed = 0;
    
//   for(ci_string::const_iterator first = s.begin(); first != s.end(); ++first) {
//     std::tr1::hash<char> hasher;
//     seed ^= hasher(std::tolower(*first)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//   }
    
//   return seed;
// }

// } // namespace boost
