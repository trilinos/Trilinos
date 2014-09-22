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

#ifndef STK_UTIL_UTIL_CI_TRAITS_H
#define STK_UTIL_UTIL_CI_TRAITS_H

#include <cctype>
#include <string>

/**
 * @brief Class <b>ignorecase_traits</b> is a character traits class that ignores case during compares.
 *
 * Replace functions of the standard char_traits<char> so that strings behave in a case-insensitive
 * way.
 *
 */
struct ignorecase_traits : public std::char_traits<char>
{
  /**
   * @brief Member function <b>eq</b> return true is c1 and c2 are equal.
   *
   * @param c1			a <b>char</b> const reference to character to compare. 
   *
   * @param c2			a <b>char</b> const reference to character to compare. 
   *
   * @return			a <b>bool</b> of true if c1 and c2 are equal
   */
  static bool eq(const char &c1, const char &c2) {
    return std::toupper(c1) == std::toupper(c2);
  }
  
  /**
   * @brief Member function <b>lt</b> return true is c1 less than c2.
   *
   * @param c1			a <b>char</b> const ...
   *
   * @param c2			a <b>char</b> const ...
   *
   * @return			a <b>bool</b> ...
   */
  static bool lt(const char &c1, const char &c2) {
    return std::toupper(c1) < std::toupper(c2);
  }
  
  /**
   * @brief Member function <b>compare</b> compares up to n characters of s1 and s2 and returns -1 if
   * s1 is less then s2, 0 if they are equal, and 1 if s1 is greater than s2.
   *
   * @param s1			a <b>char</b> const pointer to string to compare.
   *
   * @param s2			a <b>char</b> const pointer to string to compare.
   *
   * @param n			a <b>std::size_t</b> maxiumum number of character to compare. 
   *
   * @return                    an <b>int</b> value of -1 if s1 is less then s2, 0 if they are
   *                            equal, and 1 if s1 is greater than s2.
   */
  static int compare(const char *s1, const char *s2, std::size_t n);
  
  /**
   * @brief Member function <b>find</b> returns char pointer to first occurrence of character
   * <b>c</b> in first <b>n</b> characters of string <b>s</b> or 0 if not found.
   *
   * @param s			a <b>char</b> const pointer to string to search in.
   *
   * @param n			a <b>std::size_t</b> value of the maximum number of characters to
   *                            compare. 
   *
   * @param c			a <b>char</b> const reference to the character to search.
   *
   * @return			a <b>char</b> pointer to first occurrence or 0 is not found.
   */
  static const char *find(const char *s, std::size_t n, const char &c);
};

#endif // STK_UTIL_UTIL_CI_TRAITS_H
