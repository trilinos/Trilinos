// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_TransformBCNameForIOSS.hpp"
#include "Panzer_String_Utilities.hpp"
#include <algorithm>
#include <cctype>

std::string panzer_stk::transformBCNameForIOSS(std::string& name)
{
  // strip off leading and trailing whitespace just in case this comes
  // in from input file.
  panzer::trim(name);

  // replace internal whitespace with underscores and upper case with lower case.
  std::transform(name.begin(), name.end(), name.begin(),
                 [](const char c)
                 {
                   if (c == ' ')
                     return '_';
                   else
                     return char(std::tolower(c));
                 });
  return name;
}
