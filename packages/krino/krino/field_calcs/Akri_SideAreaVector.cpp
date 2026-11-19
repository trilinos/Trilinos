// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SideAreaVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

static stk::math::Vector3d compute_line_or_tri_side_area_vector(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const unsigned dim,
    const StkMeshEntities & sideNodes)
{
  STK_ThrowRequireMsg(sideNodes.size() == 2 || sideNodes.size() == 3,
      "Side area vector expects linear line or tri sides.");
  const stk::math::Vector3d x0 = get_vector_field(mesh, coordsField, sideNodes[0], dim);
  const stk::math::Vector3d x1 = get_vector_field(mesh, coordsField, sideNodes[1], dim);
  if (sideNodes.size() == 2) return crossZ(x1 - x0);
  const stk::math::Vector3d x2 = get_vector_field(mesh, coordsField, sideNodes[2], dim);
  return 0.5 * Cross(x1 - x0, x2 - x0);
}

static stk::mesh::Entity
get_element_node_not_on_side(const StkMeshEntities & elemNodes, const StkMeshEntities & sideNodes)
{
  for (auto node : elemNodes)
    if (std::find(sideNodes.begin(), sideNodes.end(), node) == sideNodes.end()) return node;
  STK_ThrowRequireMsg(false, "Element has no node off its side.");
  return stk::mesh::Entity();
}

stk::math::Vector3d
compute_line_or_tri_side_area_vector_out_of_element(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Entity side,
    const stk::mesh::Entity element)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const StkMeshEntities elemNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  STK_ThrowRequireMsg(
      elemNodes.size() == dim + 1, "Side area vector expects linear tri or tet elements.");
  const StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
  const stk::math::Vector3d areaVec =
      compute_line_or_tri_side_area_vector(mesh, coordsField, dim, sideNodes);
  const stk::math::Vector3d xOpposite =
      get_vector_field(mesh, coordsField, get_element_node_not_on_side(elemNodes, sideNodes), dim);
  const bool pointsIntoElement =
      Dot(areaVec, xOpposite - get_vector_field(mesh, coordsField, sideNodes[0], dim)) > 0.;
  return pointsIntoElement ? -areaVec : areaVec;
}

} // namespace krino
