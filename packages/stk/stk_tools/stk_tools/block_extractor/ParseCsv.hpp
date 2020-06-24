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

#ifndef PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_
#define PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace tools {

std::string strip_string(const std::string &token)
{
  std::string tmp(token);
  tmp.erase(0, tmp.find_first_not_of(" "));
  tmp.erase(tmp.find_last_not_of(" ")+1);
  return tmp;
}

std::vector<std::string> split_string(const std::string &input, const char separator)
{
  std::vector<std::string> separated;
  std::istringstream iss(input);
  std::string token;
  while(std::getline(iss, token, separator))
    separated.push_back(strip_string(token));
  return separated;
}

std::vector<std::string> get_csv(const std::string &input)
{
  return split_string(input, ',');
}

int string_to_integer(const std::string &input)
{
  int convertedInt = 0;

  try {
    convertedInt = std::stoi(input);
  }
  catch (std::exception &) {
    ThrowErrorMsg("Could not convert argument '" << input << "' to an integer");
  }

  return convertedInt;
}


std::vector<int> get_ids_from_strings(const std::vector<std::string> & input)
{
  std::vector<int> ids;
  for (const std::string & inputString : input)
  {
    const std::vector<std::string> splitInput = split_string(inputString, ':');
    ThrowRequireMsg((splitInput.size() > 0) && (splitInput.size() <= 3),
                    "Incorrect number of fields for range (" << inputString << ").  Syntax is first:last:stride");

    const int first = string_to_integer(splitInput[0]);
    const int last = (splitInput.size() > 1) ? string_to_integer(splitInput[1]) : first;
    const int stride = (splitInput.size() > 2) ? string_to_integer(splitInput[2]) : 1;

    ThrowRequireMsg((first > 0) && (last > 0) && (stride > 0) && (last >= first),
                    "Invalid integers provided for range (" << inputString << ").  Syntax is first:last:stride");

    for (int i = first; i <= last; i += stride) {
      ids.push_back(i);
    }
  }
  return ids;
}

}}

#endif /* PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_ */
