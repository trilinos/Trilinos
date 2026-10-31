// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_NodeMarker.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_Intersection_Points.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace krino {

void mark_nodes_near_edge_intersections(const stk::mesh::BulkData& mesh,
    const std::vector<IntersectionPoint> & edgeIntersections,
    FieldRef node_marker_field)
{
  const double overlap = 0.25;

  stk::mesh::field_fill(0, node_marker_field);

  auto nodeMarkerData = node_marker_field.data<int, stk::mesh::ReadWrite>();

  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);

    if (edge.crossingLocation < 0.5+overlap)
      nodeMarkerData.entity_values(edge.nodes[0])() = 1;
    if (edge.crossingLocation > 0.5-overlap)
      nodeMarkerData.entity_values(edge.nodes[1])() = 1;
  }
  stk::mesh::parallel_max(mesh, {&node_marker_field.field()});
}

bool element_has_marked_node(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem, const FieldRef nodeMarkerField)
{
  const unsigned num_nodes = mesh.num_nodes(elem);
  const stk::mesh::Entity * const elemNodes = mesh.begin_nodes(elem);
  auto nodeMarkerData = nodeMarkerField.data<int>();
  for(unsigned i=0; i < num_nodes; ++i)
  {
    const int nodeMarker = nodeMarkerData.entity_values(elemNodes[i])();
    if(nodeMarker)
      return true;
  }
  return false;
}

}
