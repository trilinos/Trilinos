//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_String_Utilities.hpp"
#include <sstream>

namespace Tempus {

void trim(std::string& str)
{
  const std::string whitespace(" \t\n");

  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos) {
    str = "";
    return;  // no content
  }

  const auto strEnd   = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  str = str.substr(strBegin, strRange);
  return;
}

void StringTokenizer(std::vector<std::string>& tokens, const std::string& str,
                     const std::string delimiters, bool trim)
{
  using std::string;

  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // grab token, trim if desired
    std::string token = str.substr(lastPos, pos - lastPos);
    if (trim) Tempus::trim(token);

    // Found a token, add it to the vector.
    tokens.push_back(token);

    if (pos == string::npos) break;

    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void TokensToDoubles(std::vector<double>& values,
                     const std::vector<std::string>& tokens)
{
  // turn tokens into doubles (its a miracle!)
  for (std::size_t i = 0; i < tokens.size(); i++) {
    double value = 0.0;
    std::stringstream ss;
    ss << tokens[i];
    ss >> value;

    values.push_back(value);
  }
}

void TokensToInts(std::vector<int>& values,
                  const std::vector<std::string>& tokens)
{
  // turn tokens into doubles (its a miracle!)
  for (std::size_t i = 0; i < tokens.size(); i++) {
    int value = 0;
    std::stringstream ss;
    ss << tokens[i];
    ss >> value;

    values.push_back(value);
  }
}
}  // namespace Tempus
