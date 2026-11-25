// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshTopologySharpFeatures.hpp>

#include <algorithm>
#include <array>
#include <map>

#include <Akri_MeshHelpers.hpp>
#include <Akri_SharpFeatureCalcs.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

namespace {

// The ordinals of the side-rank parts (other than stk's topology root parts) of a side: the
// surface it is in.
using SurfaceKey = std::vector<unsigned>;

SurfaceKey get_surface_key(const stk::mesh::BulkData & mesh, const stk::mesh::Entity side)
{
  const stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();
  SurfaceKey key;
  for (auto * part : mesh.bucket(side).supersets())
    if (part->primary_entity_rank() == sideRank && !stk::mesh::is_topology_root_part(*part))
      key.push_back(part->mesh_meta_data_ordinal());
  std::sort(key.begin(), key.end());
  return key;
}

bool side_has_selected_element(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Entity side)
{
  const stk::mesh::Entity * elems = mesh.begin_elements(side);
  for (unsigned i = 0; i < mesh.num_elements(side); ++i)
    if (elementSelector(mesh.bucket(elems[i])))
      return true;
  return false;
}

bool node_is_in_a_nodeset(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  for (auto * part : mesh.bucket(node).supersets())
    if (part->primary_entity_rank() == stk::topology::NODE_RANK && stk::io::is_part_io_part(*part) &&
        !stk::mesh::is_topology_root_part(*part))
      return true;
  return false;
}

bool edge_between_nodes_is_in_an_edge_part(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Entity otherNode)
{
  const stk::mesh::Entity * edges = mesh.begin(node, stk::topology::EDGE_RANK);
  for (unsigned i = 0; i < mesh.num_connectivity(node, stk::topology::EDGE_RANK); ++i)
  {
    const stk::mesh::Entity * edgeNodes = mesh.begin_nodes(edges[i]);
    if (edgeNodes[0] != otherNode && edgeNodes[1] != otherNode) continue;
    for (auto * part : mesh.bucket(edges[i]).supersets())
      if (part->primary_entity_rank() == stk::topology::EDGE_RANK && !stk::mesh::is_topology_root_part(*part))
        return true;
  }
  return false;
}

struct NodeSide
{
  stk::mesh::Entity side;
  SurfaceKey surface;
};

// The node's sides that are on a surface and have a selected element.
std::vector<NodeSide> get_node_surface_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Entity node)
{
  const stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();
  std::vector<NodeSide> nodeSides;
  const stk::mesh::Entity * sides = mesh.begin(node, sideRank);
  for (unsigned i = 0; i < mesh.num_connectivity(node, sideRank); ++i)
  {
    if (!side_has_selected_element(mesh, elementSelector, sides[i])) continue;
    SurfaceKey surface = get_surface_key(mesh, sides[i]);
    if (!surface.empty())
      nodeSides.push_back({sides[i], std::move(surface)});
  }
  return nodeSides;
}

stk::math::Vector3d node_location(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity node)
{
  return get_vector_field(mesh, coordsField, node, mesh.mesh_meta_data().spatial_dimension());
}

// The unit tangent of a curve through the node with these two neighbors on it, or zero if they are in
// the same direction from the node.
stk::math::Vector3d curve_tangent(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity node,
    const std::array<stk::mesh::Entity, 2> & neighbors)
{
  const stk::math::Vector3d x = node_location(mesh, coordsField, node);
  const stk::math::Vector3d tangent = (node_location(mesh, coordsField, neighbors[0]) - x).unit_vector() -
      (node_location(mesh, coordsField, neighbors[1]) - x).unit_vector();
  const double length = tangent.length();
  return (length > 0.) ? tangent / length : stk::math::Vector3d::ZERO;
}

// The area vector of a planar polygonal side from its vertices, in order.
stk::math::Vector3d side_area_vector(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity side)
{
  const unsigned numVertices = mesh.bucket(side).topology().base().num_nodes();
  const stk::mesh::Entity * sideNodes = mesh.begin_nodes(side);
  stk::math::Vector3d areaVector = stk::math::Vector3d::ZERO;
  for (unsigned i = 0; i < numVertices; ++i)
    areaVector += 0.5 * Cross(node_location(mesh, coordsField, sideNodes[i]), node_location(mesh, coordsField, sideNodes[(i + 1) % numVertices]));
  return areaVector;
}

SmoothingDirection fixed_smoothing_direction(const stk::mesh::Entity node)
{
  return {node, SmoothingDirection::FIXED, stk::math::Vector3d::ZERO};
}

// The point off the edge of a side: its vertex off the edge for a triangle, else the average of its
// vertices off the edge.
stk::math::Vector3d side_point_off_edge(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity side,
    const stk::mesh::Entity edgeNode0, const stk::mesh::Entity edgeNode1)
{
  const unsigned numVertices = mesh.bucket(side).topology().base().num_nodes();
  const stk::mesh::Entity * sideNodes = mesh.begin_nodes(side);
  stk::math::Vector3d sum = stk::math::Vector3d::ZERO;
  unsigned numOffEdge = 0;
  for (unsigned i = 0; i < numVertices; ++i)
  {
    if (sideNodes[i] == edgeNode0 || sideNodes[i] == edgeNode1) continue;
    sum += node_location(mesh, coordsField, sideNodes[i]);
    ++numOffEdge;
  }
  return sum / numOffEdge;
}

bool sides_meet_at_sharp_angle_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double cosFeatureAngle,
    const stk::mesh::Entity edgeNode0, const stk::mesh::Entity edgeNode1, const std::array<stk::mesh::Entity, 2> & sides)
{
  return is_angle_sharper_than_feature_angle(
      cosine_of_angle_between_triangles_at_edge(node_location(mesh, coordsField, edgeNode0),
          node_location(mesh, coordsField, edgeNode1),
          side_point_off_edge(mesh, coordsField, sides[0], edgeNode0, edgeNode1),
          side_point_off_edge(mesh, coordsField, sides[1], edgeNode0, edgeNode1)),
      cosFeatureAngle);
}

SmoothingDirection smoothing_direction_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity node,
    const std::vector<NodeSide> & nodeSides, const std::optional<double> cosFeatureAngle)
{
  // The sides using each edge from the node, by the edge's other node.
  std::map<stk::mesh::Entity, std::vector<const NodeSide *>> edgeSides;
  for (auto && nodeSide : nodeSides)
  {
    const stk::topology sideTopology = mesh.bucket(nodeSide.side).topology();
    const stk::mesh::Entity * sideNodes = mesh.begin_nodes(nodeSide.side);
    for (unsigned iEdge = 0; iEdge < sideTopology.num_edges(); ++iEdge)
    {
      std::array<unsigned, 3> edgeOrdinals;
      sideTopology.edge_node_ordinals(iEdge, edgeOrdinals.data());
      const stk::mesh::Entity edgeNode0 = sideNodes[edgeOrdinals[0]];
      const stk::mesh::Entity edgeNode1 = sideNodes[edgeOrdinals[1]];
      if (edgeNode0 == node) edgeSides[edgeNode1].push_back(&nodeSide);
      else if (edgeNode1 == node) edgeSides[edgeNode0].push_back(&nodeSide);
    }
  }

  std::vector<stk::mesh::Entity> featureEdgeNeighbors;
  for (auto && [neighbor, sides] : edgeSides)
    if (sides.size() != 2 || sides[0]->surface != sides[1]->surface ||
        edge_between_nodes_is_in_an_edge_part(mesh, node, neighbor) ||
        (cosFeatureAngle && sides_meet_at_sharp_angle_3D(mesh, coordsField, *cosFeatureAngle, node, neighbor, {{sides[0]->side, sides[1]->side}})))
      featureEdgeNeighbors.push_back(neighbor);

  if (featureEdgeNeighbors.empty())
  {
    // Every side of the node is in one surface around it (unless separate fans of a surface touch at it).
    for (auto && nodeSide : nodeSides)
      if (nodeSide.surface != nodeSides[0].surface) return fixed_smoothing_direction(node);
    stk::math::Vector3d normal = stk::math::Vector3d::ZERO;
    for (auto && nodeSide : nodeSides)
    {
      const stk::math::Vector3d areaVector = side_area_vector(mesh, coordsField, nodeSide.side);
      normal += (Dot(areaVector, normal) < 0.) ? -areaVector : areaVector;
    }
    const double length = normal.length();
    if (length == 0.) return fixed_smoothing_direction(node);
    return {node, SmoothingDirection::SMOOTH_IN_PLANE, normal / length};
  }

  if (featureEdgeNeighbors.size() == 2)
  {
    const stk::math::Vector3d tangent = curve_tangent(mesh, coordsField, node, {{featureEdgeNeighbors[0], featureEdgeNeighbors[1]}});
    if (tangent.length() == 0.) return fixed_smoothing_direction(node);
    return {node, SmoothingDirection::SMOOTH_ALONG_TANGENT, tangent};
  }

  return fixed_smoothing_direction(node);
}

SmoothingDirection smoothing_direction_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity node,
    const std::vector<NodeSide> & nodeSides, const std::optional<double> cosFeatureAngle)
{
  if (nodeSides.size() != 2 || nodeSides[0].surface != nodeSides[1].surface) return fixed_smoothing_direction(node);

  std::array<stk::mesh::Entity, 2> neighbors;
  for (unsigned i = 0; i < 2; ++i)
  {
    const stk::mesh::Entity * segmentNodes = mesh.begin_nodes(nodeSides[i].side);
    neighbors[i] = (segmentNodes[0] == node) ? segmentNodes[1] : segmentNodes[0];
  }
  if (cosFeatureAngle &&
      is_angle_sharper_than_feature_angle(cosine_of_angle_between_segments_at_vertex(node_location(mesh, coordsField, node),
          node_location(mesh, coordsField, neighbors[0]), node_location(mesh, coordsField, neighbors[1])), *cosFeatureAngle))
    return fixed_smoothing_direction(node);
  const stk::math::Vector3d tangent = curve_tangent(mesh, coordsField, node, neighbors);
  if (tangent.length() == 0.) return fixed_smoothing_direction(node);
  return {node, SmoothingDirection::SMOOTH_ALONG_TANGENT, tangent};
}

}

std::vector<SmoothingDirection> compute_smoothing_directions_from_mesh_topology(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const std::optional<double> cosFeatureAngle)
{
  STK_ThrowRequireMsg(mesh.parallel_size() == 1 || mesh.is_automatic_aura_on(),
      "Finding the smoothing directions on more than one processor needs the aura, for the sides of the shared nodes.");

  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const stk::mesh::Selector ownedOrSharedSelectedNodes = elementSelector & (meta.locally_owned_part() | meta.globally_shared_part());

  std::vector<SmoothingDirection> smoothingDirections;
  for (auto * bucket : mesh.get_buckets(stk::topology::NODE_RANK, ownedOrSharedSelectedNodes))
  {
    for (auto node : *bucket)
    {
      const std::vector<NodeSide> nodeSides = get_node_surface_sides(mesh, elementSelector, node);
      const bool isInANodeset = node_is_in_a_nodeset(mesh, node);
      if (nodeSides.empty() && !isInANodeset) continue;

      if (isInANodeset)
        smoothingDirections.push_back(fixed_smoothing_direction(node));
      else if (meta.spatial_dimension() == 3)
        smoothingDirections.push_back(smoothing_direction_3D(mesh, coordsField, node, nodeSides, cosFeatureAngle));
      else
        smoothingDirections.push_back(smoothing_direction_2D(mesh, coordsField, node, nodeSides, cosFeatureAngle));
    }
  }

  std::sort(smoothingDirections.begin(), smoothingDirections.end(),
      [&mesh](const SmoothingDirection & a, const SmoothingDirection & b) { return mesh.identifier(a.node) < mesh.identifier(b.node); });
  return smoothingDirections;
}

}
