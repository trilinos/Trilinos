// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_INTERFACE_GEOMETRY_UNIT_TESTS_AKRI_UNIT_SHARPFEATUREINTERSECTIONPOINTCHECKS_HPP_
#define KRINO_KRINO_INTERFACE_GEOMETRY_UNIT_TESTS_AKRI_UNIT_SHARPFEATUREINTERSECTIONPOINTCHECKS_HPP_

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_ElementCutterUtils.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_ParametricElementLibrary.hpp>
#include <Akri_SharpFeatureIntersectionPoints.hpp>

namespace krino {

// Checks of the sharp-feature points of the elements of the mesh of a fixture, by element and node
// index: BASE is a mesh fixture with mMesh, get_assigned_element_for_index and
// get_assigned_node_for_index.
template <typename BASE>
class SharpFeaturePointChecks : public BASE
{
public:
  // The points of one sharp-feature rank in an element, as the cutting of that rank gets them: its
  // interior and face points, reduced by tol, in canonical order, without repeats.
  std::vector<IntersectionPoint> element_points(const InterfaceGeometry & geometry,
      const size_t elementIndex,
      const int codimension,
      const double tol)
  {
    const stk::mesh::Entity element = this->get_assigned_element_for_index(elementIndex);
    ParametricElementLibrary & parametricElemLibrary =
        AuxMetaData::get(this->mMesh.mesh_meta_data()).parametric_element_library();
    const ParametricElement & parametricElement =
        parametricElemLibrary.get_parametric_element(this->mMesh.bucket(element).topology());
    const std::vector<stk::mesh::Entity> elementNodes(
        this->mMesh.begin_nodes(element), this->mMesh.end_nodes(element));

    geometry.prepare_to_decompose_elements(this->mMesh, NodeToCapturedDomainsMap());
    const std::unique_ptr<ElementCutter> cutter = geometry.build_element_cutter(
        this->mMesh, element, [](const std::array<unsigned, 4> &) { return true; });

    std::vector<IntersectionPoint> points;
    append_intersection_points_from_element_interior(
        parametricElement, elementNodes, *cutter, keep_all_intersection_points_filter(), points);
    for (int iFace = 0; iFace < 4; ++iFace)
      append_intersection_points_from_element_face(parametricElemLibrary,
          parametricElement,
          elementNodes,
          iFace,
          *cutter,
          keep_all_intersection_points_filter(),
          points);
    return reduced_intersection_points_of_codimension(this->mMesh, points, codimension, tol);
  }

  // The points whose stencil is the nodes with these indices.
  std::vector<const IntersectionPoint *>
  points_on(const std::vector<IntersectionPoint> & points, const std::vector<size_t> & nodeIndices)
  {
    std::vector<stk::mesh::Entity> nodes = get_nodes(nodeIndices);
    std::sort(nodes.begin(), nodes.end());
    std::vector<const IntersectionPoint *> pointsOn;
    for (auto & point : points)
    {
      std::vector<stk::mesh::Entity> pointNodes = point.get_nodes();
      std::sort(pointNodes.begin(), pointNodes.end());
      if (pointNodes == nodes) pointsOn.push_back(&point);
    }
    return pointsOn;
  }

  // Each of the elements with these indices has one point with this stencil, the same in each.
  // Returns the weight of the first node.
  double expect_same_point_on_in_elements(const InterfaceGeometry & geometry,
      const std::vector<size_t> & elementIndices,
      const std::vector<size_t> & nodeIndices,
      const int codimension,
      const double tol)
  {
    std::vector<IntersectionPoint> pointsInElements;
    for (size_t elementIndex : elementIndices)
    {
      const std::vector<IntersectionPoint> points =
          element_points(geometry, elementIndex, codimension, tol);
      const std::vector<const IntersectionPoint *> pointsOn = points_on(points, nodeIndices);
      EXPECT_EQ(1u, pointsOn.size()) << "element " << elementIndex;
      if (pointsOn.size() == 1) pointsInElements.push_back(*pointsOn[0]);
    }
    if (pointsInElements.size() != elementIndices.size()) return 0.;
    for (auto & point : pointsInElements)
    {
      EXPECT_TRUE(is_same_intersection_point(pointsInElements[0], point, 1.e-12));
    }
    return weight_of_node(pointsInElements[0], nodeIndices.front());
  }

  double weight_of_node(const IntersectionPoint & point, const size_t nodeIndex) const
  {
    for (size_t i = 0; i < point.get_nodes().size(); ++i)
      if (point.get_nodes()[i] == this->get_assigned_node_for_index(nodeIndex))
        return point.get_weights()[i];
    return 0.;
  }

private:
  std::vector<stk::mesh::Entity> get_nodes(const std::vector<size_t> & nodeIndices) const
  {
    std::vector<stk::mesh::Entity> nodes;
    for (size_t nodeIndex : nodeIndices)
      nodes.push_back(this->get_assigned_node_for_index(nodeIndex));
    return nodes;
  }
};

} // namespace krino

#endif /* KRINO_KRINO_INTERFACE_GEOMETRY_UNIT_TESTS_AKRI_UNIT_SHARPFEATUREINTERSECTIONPOINTCHECKS_HPP_ */
