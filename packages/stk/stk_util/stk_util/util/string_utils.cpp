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

#include "stk_util/util/string_utils.hpp"
#include <cstddef>  // for size_t
#include <string>   // for string, operator+, allocator, char_traits
#include <cctype>
#include <algorithm>
#include <sstream>

namespace stk {

std::string angle_it(const std::string &s)
{
    return "<" + s + ">";
}
std::string bracket_it(const std::string &s)
{
    return "[" + s + "]";
}
std::string dash_it(const std::string &s)
{
    return "--" + s;
}

std::string rm_dashes(const std::string &s)
{
  size_t pos = s.find_first_not_of(" -");
  if (pos != std::string::npos) {
    return s.substr(pos);
  }
  return s;
}

std::string get_substring_before_comma(const std::string& s)
{
  size_t pos = s.find(",");
  return s.substr(0,pos);
}

std::string get_substring_after_comma(const std::string& s)
{
  size_t pos = s.find(",");
  if (pos != std::string::npos) {
    return s.substr(pos+1);
  }
  return "";
}

std::string tailname(const std::string &filename)
{
  size_t ind = filename.find_last_of("/", filename.size());
  if (ind != std::string::npos) {
    return filename.substr(ind+1, filename.size());
  }
  return filename; // No path, just return the filename
}

std::string basename(const std::string &filename)
{
  std::string tail = tailname(filename);

  // Strip off the extension
  size_t ind = tail.find_last_of('.', tail.size());
  if (ind != std::string::npos) {
    return tail.substr(0,ind);
  }
  return tail;
}

std::vector<std::string> make_vector_of_strings(const std::string& inputString, char separator, int maxStringLength)
{
  std::vector<std::string> vecStr;
  int processedLen = 0;
  int descrSize = inputString.size();
  while(processedLen < descrSize) {
    size_t nextLineBreak = inputString.find('\n', processedLen);
    if (nextLineBreak != std::string::npos && static_cast<int>(nextLineBreak) <= processedLen+maxStringLength) {
      vecStr.push_back(inputString.substr(processedLen, nextLineBreak-processedLen));
      processedLen = nextLineBreak+1;
      continue;
    }

    int remaining = inputString.size() - processedLen;
    if (remaining <= maxStringLength) {
//std::cerr<<"remaining: "<<remaining<<", maxStringLength: "<<maxStringLength<<std::endl;
      vecStr.push_back(inputString.substr(processedLen));
      processedLen += remaining;
    }
    else {
      size_t pos = inputString.rfind(separator, processedLen + maxStringLength);
      if (pos == std::string::npos || pos < static_cast<unsigned>(processedLen)) {
        pos = processedLen + maxStringLength;
      }
      int numChars = pos - processedLen;
//std::cerr<<"pos: "<<pos<<", processedLen: "<<processedLen<<", substr: '"<<inputString.substr(processedLen, numChars)<<"'"<<std::endl;
      vecStr.push_back(inputString.substr(processedLen, numChars));
      if (inputString[pos] == separator) {
        pos = inputString.find_first_not_of(separator, pos);
        if (pos == std::string::npos) {
          break;
        }
      }
      processedLen = pos;
    }
  }
  return vecStr;
}

bool string_starts_with(const std::string & queryString, const std::string & prefix)
{
  return (queryString.rfind(prefix, 0) == 0);
}

std::string ltrim_string(std::string s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  [](unsigned char ch) { return !std::isspace(ch); }
                                  ));
  return s;
}

std::string rtrim_string(std::string s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); }
                       ).base(), s.end());
  return s;
}

std::string trim_string(std::string s) {
  return ltrim_string(rtrim_string(s));
}

std::vector<std::string> split_string(const std::string & input, const char separator)
{
  std::vector<std::string> separated;
  std::istringstream iss(input);
  std::string token;
  while (std::getline(iss, token, separator)) {
    separated.push_back(token);
  }
  if (input.back() == separator) {
    separated.push_back("");
  }
  return separated;
}

std::vector<std::string> split_csv_string(const std::string & input)
{
  return split_string(input, ',');
}

} // namespace stk

