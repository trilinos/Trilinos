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

std::string panzer_stk::transformBCNameForIOSS(std::string& name, const bool make_lower_case)
{
  // strip off leading and trailing whitespace just in case this comes
  // in from input file.
  panzer::trim(name);

  std::transform(name.begin(), name.end(), name.begin(),
                 [&make_lower_case](const char c)
                 {
                   if (c == ' ')
                     return '_';
                   else if (make_lower_case)
                     return char(std::tolower(c));

                   return char(c);
                 });
  return name;
}
