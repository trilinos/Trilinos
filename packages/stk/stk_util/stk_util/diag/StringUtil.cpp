/*
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
 */

#include "stk_util/diag/StringUtil.hpp"

#include <cstdlib>  // for free
#include <iomanip>  // for operator<<, setprecision
#include <sstream>  // for operator<<, basic_ostream, ostream, basic_ostream::operator<<, ostri...
#include <string>   // for basic_string, string, basic_string<>::const_iterator, char_traits

#if __GNUC__ >= 3
#include <cxxabi.h>  // for __cxa_demangle
#endif

//----------------------------------------------------------------------

namespace sierra {

int
case_strcmp(
  const char *c1,
  const char *c2)
{
  for ( ; ; c1++, c2++) {
    if ( std::tolower(*c1) != std::tolower(*c2) )
      return ( std::tolower(*c1) - std::tolower(*c2) ) ;
    if (*c1 == '\0')
      return 0 ;
  }
}

const char *
case_strstr(
  const char *s,
  const char *find)
{
  if (!*find)
    return s;

  const char *cp = s;
  while (*cp) {
    const char *t1 = cp;
    const char *t2 = find;

    for ( ; std::tolower(*t1) && std::tolower(*t2) && !(std::tolower(*t1) - std::tolower(*t2)); ++t1, ++t2)
      ;

    if (!*t2)
      return cp;

    cp++;
  }

  return nullptr;
}


std::string
title(
  const std::string &s)
{
  std::string t(s);

  bool all_upper = true;
  bool all_lower = true;

  bool next_upper = true;
  for (std::string::iterator c = t.begin(); c != t.end(); ++c) {
    all_upper &= (*c == std::toupper(*c));
    all_lower &= (*c == std::tolower(*c));
    if (next_upper)
      *c = std::toupper(*c);
    else
      *c = std::tolower(*c);
    next_upper = !isalpha(*c);
  }
  
  if (all_upper || all_lower) 
    return t;
  else
    return s;
}


template std::string to_string<double>(const double &);
template std::string to_string<float>(const float &);
template std::string to_string<int>(const int &);
template std::string to_string<unsigned>(const unsigned &);
template std::string to_string<long>(const long &);
template std::string to_string<unsigned long>(const unsigned long &);

std::string
to_string(
  const double &r,
  int precision)
{
  std::ostringstream os;
  os << std::setprecision(precision) << r;
  return std::string(os.str());
}


std::ostream &
object_phrase::print(
  std::ostream &os) const
{
  if (m_n == 0)
    os << m_plural << " no " << m_noun << "s";
  else if (m_n == 1)
    os << m_singular << " 1 " << m_noun;
  else
    os << m_plural << " " << m_n << " " << m_noun << "s";

  return os;
}

namespace {

std::string::const_iterator
find_next_char(
  std::string::const_iterator p,
  std::string::const_iterator end,
  char c)
{
  while (p != end && *p != c)
    p++;
  return p;
}

std::string::const_iterator
find_next_not_char(
  std::string::const_iterator p,
  std::string::const_iterator end,
  char c)
{
  while (p != end && *p == c)
    p++;
  return p;
}

inline std::string::const_iterator find_next_space(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, ' ');
}

inline std::string::const_iterator find_next_endl(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, '\n');
}

inline std::string::const_iterator find_next_nonspace(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_not_char(p, end, ' ');
}

auto find_next_space_after_nonspace(std::string::const_iterator p, std::string::const_iterator end)
{
  return find_next_space(find_next_nonspace(p, end), end);
}

} // namespace <null>

std::string word_wrap(
    const std::string &s, unsigned int max_line_length, const std::string &prefix, const std::string &prefix_first_line)
{
  std::string t;
  const std::string *u = &prefix_first_line;

  std::string::const_iterator line_start, next_space, line_end, next_newline;
  line_start = line_end = s.begin();

  const auto line_length = [&line_start, &u](std::string::const_iterator _line_end) {
    return std::distance(line_start, _line_end) + u->size();
  };

  while (line_end != s.end()) {
    line_end = find_next_space_after_nonspace(line_start, s.end());
    do {
      next_space = find_next_space_after_nonspace(line_end, s.end());
      if (line_length(next_space) > max_line_length) {
        break;
      }
      line_end = next_space;
    } while (line_end != s.end());

    next_newline = find_next_endl(line_start, s.end());
    if ((next_newline < line_end) || (line_length(next_newline) < max_line_length)) {
      line_end = next_newline;
    }

    t.append(*u).append(line_start, line_end).append("\n");

    if (line_end == next_newline)
      u = &prefix_first_line;
    else
      u = &prefix;

    line_start = line_end + 1;
  }

  return t;
}

#ifdef SIERRA_USE_PLATFORM_DEMANGLER

  #if defined(__GNUC__)

    #if (__GNUC__ == 3)
      std::string
      demangle(const char * symbol)
     {
       std::string s;
       int status;

       char *demangled_symbol = abi::__cxa_demangle(symbol, 0, 0, &status);

       if (demangled_symbol) {
         s = std::string(demangled_symbol);
         free(demangled_symbol);
       }

       if (status != 0)
         s = std::string(symbol);

       return s;
     }

    #elif (__GNUC__ == 4)
      std::string
      demangle(const char * symbol)
      {
        std::string s;

        int status = -1;

        char *demangled_symbol = __cxxabiv1::__cxa_demangle(symbol, 0, 0, &status);

        if (demangled_symbol) {
          s = std::string(demangled_symbol);
          free(demangled_symbol);
        }

        if (status != 0)
          s = std::string(symbol);

        return s;
      }

    #elif (__GNUC__ >= 5)
      std::string
      demangle(const char * symbol)
      {
        std::string s;

        int status = -1;

        char *demangled_symbol = abi::__cxa_demangle(symbol, 0, 0, &status);

        if (demangled_symbol) {
          s = std::string(demangled_symbol);
          free(demangled_symbol);
        }

        if (status != 0)
          s = std::string(symbol);

        return s;
      }
      
    #endif // (__GNUC__ == 3)

#else

std::string demangle(const char *symbol) {
  return symbol;
}

#endif // defined(__GNUC__)

#else

const char *demangle(const char *symbol) {
  return symbol;
}

#endif // SIERRA_USE_PLATFORM_DEMANGLER

} // namespace sierra
