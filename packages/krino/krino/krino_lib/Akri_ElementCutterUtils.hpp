// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_
#define KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_
#include <stk_math/StkVector.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <vector>

namespace krino {

class ElementCutter;
class InterfaceGeometry;

void append_intersection_points_from_element_interior(const MasterElement & masterElement,
    const std::vector<stk::mesh::Entity> & nodes,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints);

void append_intersection_points_from_within_element_and_owned_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const stk::mesh::Entity element,
    const ElementCutter & elementCutter,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints);

void append_closest_point_intersections(std::vector<IntersectionPoint> & intersectionPoints,
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & parentElementSelector,
    const FieldRef coordsField,
    const IntersectionPointFilter & intersectionPointFilter);

template <typename Real>
std::array<int, 4> sort_quad_vertices(const std::vector<stk::math::Vec<Real, 3>> & nodeCoords)
{
  STK_ThrowRequire(nodeCoords.size() == 4);
  std::array<int, 4> sortedOrder{0, 1, 2, 3};

  const stk::math::Vec<Real, 3> n0n1 = nodeCoords[1] - nodeCoords[0];
  const stk::math::Vec<Real, 3> n0n2 = nodeCoords[2] - nodeCoords[0];
  const stk::math::Vec<Real, 3> n0n3 = nodeCoords[3] - nodeCoords[0];

  if(stk::math::Dot(stk::math::Cross(n0n1, n0n2), stk::math::Cross(n0n1, n0n3)) < 0)
  {
    sortedOrder[2] = 1;
    sortedOrder[1] = 2;
  }
  else if(stk::math::Dot(stk::math::Cross(n0n3, n0n1), stk::math::Cross(n0n3, n0n2)) < 0)
  {
    sortedOrder[2] = 3;
    sortedOrder[3] = 2;
  }
  return sortedOrder;
}

template <typename Real>
std::vector<std::array<stk::math::Vec<Real, 3>, 3>> make_tris_from_intersections(const std::vector<stk::math::Vec<Real, 3>> & nodeCoords)
{
  STK_ThrowRequire(nodeCoords.size() == 3 || nodeCoords.size() == 4);
  std::vector<std::array<stk::math::Vec<Real, 3>, 3>> tris;
  if(nodeCoords.size() == 3) tris.push_back({nodeCoords[0], nodeCoords[1], nodeCoords[2]});
  else
  {
    const auto order = sort_quad_vertices(nodeCoords);
    const Real l0 = nodeCoords[order[0]].length_squared();
    const Real l1 = nodeCoords[order[1]].length_squared();
    const Real l2 = nodeCoords[order[2]].length_squared();
    const Real l3 = nodeCoords[order[3]].length_squared();
    const Real zeroTwoNorm = l0 < l2 ? l0 : l2;
    const Real oneThreeNorm = l1 < l3 ? l1 : l3;
    if(zeroTwoNorm <= oneThreeNorm)
    {
      tris.push_back({nodeCoords[order[0]], nodeCoords[order[1]], nodeCoords[order[2]]});
      tris.push_back({nodeCoords[order[0]], nodeCoords[order[3]], nodeCoords[order[2]]});
    }
    else 
    {
      tris.push_back({nodeCoords[order[0]], nodeCoords[order[1]], nodeCoords[order[3]]});
      tris.push_back({nodeCoords[order[1]], nodeCoords[order[2]], nodeCoords[order[3]]});
    }
  }
  return tris;
}

}

#endif /* KRINO_INCLUDE_AKRI_ELEMENTCUTTERUTILS_H_ */
