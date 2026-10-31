// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Utility.hpp>

#include <algorithm>
#include <string>

namespace krino {
namespace utility {

bool has_wildcard(const std::string & name)
{
  return name.find('*') != std::string::npos;
}

bool is_supported_wildcard_pattern(const std::string & pattern)
{
  return pattern.find('?') == std::string::npos &&
      std::count(pattern.begin(), pattern.end(), '*') <= 1;
}

bool is_wildcard_match(const std::string & pattern, const std::string & name)
{
  const std::size_t starPos = pattern.find('*');

  if (starPos == std::string::npos)
    return pattern == name;

  const std::string prefix = pattern.substr(0, starPos);
  const std::string suffix = pattern.substr(starPos + 1);

  if (name.size() < prefix.size() + suffix.size())
    return false;

  if (!prefix.empty() && name.compare(0, prefix.size(), prefix) != 0)
    return false;

  if (!suffix.empty() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
    return false;

  return true;
}

} // namespace utility
} // namespace krino
