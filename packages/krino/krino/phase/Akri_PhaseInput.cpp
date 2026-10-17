// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_PhaseInput.hpp>

#include <algorithm>
#include <cctype>

namespace krino {

Surface_Identifier SurfaceNames::get_identifier(const std::string & name)
{
  std::string upperName = name;
  std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
  for (unsigned i=0; i<myNames.size(); ++i)
    if (upperName == myNames[i])
      return Surface_Identifier(i);

  myNames.push_back(upperName);
  return Surface_Identifier(myNames.size()-1);
}

} // namespace krino
