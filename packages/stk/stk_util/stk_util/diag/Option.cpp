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

#include <sstream>
#include <iomanip>
#include <map>

#include <stk_util/diag/Option.hpp>
#include <iostream>

namespace stk {
namespace diag {

OptionMaskParser::Mask
OptionMaskParser::parse(
  const char *          mask) const
{
  if (mask) {
    const std::string mask_string(mask);

    m_status = true;

    std::string::const_iterator it0 = mask_string.begin();
    std::string::const_iterator it1;
    std::string::const_iterator it2;
    std::string::const_iterator it3;
    do {
      // Trim preceeding spaces
      while (it0 != mask_string.end() && *it0 == ' ')
        it0++;

      if (it0 == mask_string.end())
        break;

      for (it1 = it0; it1 != mask_string.end(); ++it1) {
        if (*it1 == '(' || *it1 == ':' || *it1 == ',')
          break;
      }

      // Trim trailing spaces
      it2 = it1;
      while (it2 != it0 && *(it2 - 1) == ' ')
        --it2;

      std::string name(it0, it2);

      // Get argument list
      if (*it1 == '(') {
        it2 = it1 + 1;

        // Trim preceeding spaces
        while (it2 != mask_string.end() && *it2 == ' ')
          ++it2;

        int paren_count = 0;

        for (; it1 != mask_string.end(); ++it1) {
          if (*it1 == '(')
            ++paren_count;
          else if (*it1 == ')') {
            --paren_count;
            if (paren_count == 0)
              break;
          }
        }
        it3 = it1;

        // Trim trailing spaces
        while (it3 != it2 && *(it3 - 1) == ' ')
          --it3;

        // Find next argument start
        for (; it1 != mask_string.end(); ++it1)
          if (*it1 == ':' || *it1 == ',')
            break;
      }
      else
        it2 = it3 = it1;

      const std::string arg(it2, it3);

      parseArg(name, arg);

      it0 = it1 + 1;
    } while (it1 != mask_string.end());
  }

  return m_optionMask;
}


void
OptionMaskParser::parseArg(
  const std::string &  name,
  const std::string &  arg) const
{
  OptionMaskNameMap::const_iterator mask_entry = m_optionMaskNameMap.find(name);

  if (mask_entry != m_optionMaskNameMap.end()) m_optionMask |= (*mask_entry).second.m_mask;
  else {
    Mask  mask_hex = 0;
    std::istringstream mask_hex_stream(name.c_str());
    if (mask_hex_stream >> std::resetiosflags(std::ios::basefield) >> mask_hex)
      m_optionMask |= mask_hex;
    else
      m_status = false;
  }
}


std::ostream &
OptionMaskParser::describe(
  std::ostream &    os) const
{
  os << "Specify a comma separated list of:" << std::endl;
  for (OptionMaskNameMap::const_iterator it = m_optionMaskNameMap.begin(); it != m_optionMaskNameMap.end(); ++it)
    (*it).second.describe(os);

  return os;
}


std::ostream &
OptionMaskName::describe(
  std::ostream &  os) const
{
  return os << "  " << std::left << std::setw(20) << m_name << "\t" << m_description << std::endl;
}

} // namespace diag
} // namespace stk
