// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ELEMENT_INTERSECTIONS_H_
#define KRINO_INCLUDE_AKRI_ELEMENT_INTERSECTIONS_H_

#include <stk_math/StkVector.hpp>
#include <stk_topology/topology.hpp>
#include <vector>

namespace krino {

class Element_Cutter;

struct ElementIntersection
{
  ElementIntersection(const stk::math::Vector3d & coords, const std::vector<int> & domains)
  : parametricCoords(coords),
    sortedDomains(domains) {}

  stk::math::Vector3d parametricCoords;
  std::vector<int> sortedDomains;
};

std::ostream & operator<<(std::ostream & os, const ElementIntersection& elementIntersection);

}

#endif /* KRINO_INCLUDE_AKRI_ELEMENT_INTERSECTIONS_H_ */
