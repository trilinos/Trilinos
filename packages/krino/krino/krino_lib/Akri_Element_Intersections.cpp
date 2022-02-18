// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Cutting_Surface.hpp>
#include <Akri_Element_Cutter.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_Plane_Intersections.hpp>
#include <stk_util/util/SortAndUnique.hpp>


namespace krino {

std::ostream & operator<<(std::ostream & os, const ElementIntersection & elementIntersection)
{
  os << elementIntersection.parametricCoords << ", interfaces={ ";
  for (int domain : elementIntersection.sortedDomains)
    os << domain << " ";
  os << "}";
  return os;
}

}

