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
#include <Akri_NodalMotionConstraint.hpp>
#include <Akri_SideAreaVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino
{

static stk::mesh::Entity get_selected_element_of_side(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity side,
    const stk::mesh::Selector & elementSelector)
{
  for (auto elem : StkMeshEntities{mesh.begin_elements(side), mesh.end_elements(side)})
    if (elementSelector(mesh.bucket(elem))) return elem;
  return stk::mesh::Entity();
}

namespace
{

struct BoundarySideAtNode
{
  stk::mesh::Entity side;
  stk::math::Vector3d areaVector;
  std::vector<stk::mesh::Entity> otherNodes;
};

struct Contact
{
  std::vector<size_t> walls;
};

size_t find_root(std::vector<size_t> & parent, size_t i)
{
  while (parent[i] != i)
  {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }
  return i;
}

bool have_common_node(
    const std::vector<stk::mesh::Entity> & a, const std::vector<stk::mesh::Entity> & b)
{
  for (auto node : a)
    if (std::find(b.begin(), b.end(), node) != b.end()) return true;
  return false;
}

bool contains(const std::vector<size_t> & walls, const size_t wall)
{
  return std::find(walls.begin(), walls.end(), wall) != walls.end();
}

} // namespace

NodalMotionConstraint
compute_zero_levelset_boundary_node_constraint(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef levelSetField,
    const stk::mesh::Entity node,
    const stk::mesh::Selector & boundarySideSelector,
    const stk::mesh::Selector & elementSelector,
    const std::function<bool(stk::mesh::Entity side0, stk::mesh::Entity side1)> &
        areSidesInSameWall)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();

  std::vector<BoundarySideAtNode> sides;
  for (auto side : StkMeshEntities{mesh.begin(node, sideRank), mesh.end(node, sideRank)})
  {
    if (!boundarySideSelector(mesh.bucket(side))) continue;
    const stk::mesh::Entity elem = get_selected_element_of_side(mesh, side, elementSelector);
    if (!mesh.is_valid(elem)) continue;
    BoundarySideAtNode s{side,
        compute_line_or_tri_side_area_vector_out_of_element(mesh, coordsField, side, elem),
        {}};
    for (auto sideNode : StkMeshEntities{mesh.begin_nodes(side), mesh.end_nodes(side)})
      if (sideNode != node) s.otherNodes.push_back(sideNode);
    sides.push_back(s);
  }
  if (sides.empty()) return {};

  // Walls: sides sharing an edge from the node (in 2D, the node) unless areSidesInSameWall says
  // not.
  std::vector<size_t> parent(sides.size());
  for (size_t i = 0; i < sides.size(); ++i)
    parent[i] = i;
  for (size_t i = 0; i < sides.size(); ++i)
    for (size_t j = i + 1; j < sides.size(); ++j)
      if ((dim == 2 || have_common_node(sides[i].otherNodes, sides[j].otherNodes)) &&
          areSidesInSameWall(sides[i].side, sides[j].side))
        parent[find_root(parent, i)] = find_root(parent, j);
  std::vector<size_t> wallOfSide(sides.size());
  std::vector<size_t> rootToWall(sides.size(), sides.size());
  std::vector<stk::math::Vector3d> wallAreaVector;
  for (size_t i = 0; i < sides.size(); ++i)
  {
    const size_t root = find_root(parent, i);
    if (rootToWall[root] == sides.size())
    {
      rootToWall[root] = wallAreaVector.size();
      wallAreaVector.push_back(stk::math::Vector3d::ZERO);
    }
    wallOfSide[i] = rootToWall[root];
    wallAreaVector[wallOfSide[i]] += sides[i].areaVector;
  }

  auto lsData = levelSetField.data<double>();
  auto levelSet = [&](const stk::mesh::Entity n) { return lsData.entity_values(n)(); };
  std::vector<bool> wallHasNegativeNeighbor(wallAreaVector.size(), false);
  for (size_t i = 0; i < sides.size(); ++i)
    for (auto other : sides[i].otherNodes)
      if (levelSet(other) < 0.) wallHasNegativeNeighbor[wallOfSide[i]] = true;

  auto surface_in_wall = [&](const size_t wall) -> NodalMotionConstraint
  {
    if (wallAreaVector[wall].length_squared() == 0.) return {};
    return {NodalMotionConstraint::SURFACE, wallAreaVector[wall].unit_vector()};
  };

  if (dim == 2)
  {
    std::vector<size_t> negativeWalls;
    for (size_t w = 0; w < wallAreaVector.size(); ++w)
      if (wallHasNegativeNeighbor[w]) negativeWalls.push_back(w);
    if (negativeWalls.size() == 1) return surface_in_wall(negativeWalls[0]);
    if (wallAreaVector.size() == 1) return surface_in_wall(0);
    return {};
  }

  // Contacts: toward the cut on a side's opposite edge, or along a side edge to another node with
  // zero level set (in each wall containing that edge).
  std::vector<Contact> contacts;
  std::vector<stk::mesh::Entity> zeroNeighbors;
  for (size_t i = 0; i < sides.size(); ++i)
  {
    const auto & others = sides[i].otherNodes;
    if (others.size() == 2 && levelSet(others[0]) * levelSet(others[1]) < 0.)
      contacts.push_back({{wallOfSide[i]}});
    for (auto other : others)
      if (levelSet(other) == 0. &&
          std::find(zeroNeighbors.begin(), zeroNeighbors.end(), other) == zeroNeighbors.end())
        zeroNeighbors.push_back(other);
  }
  for (auto zeroNeighbor : zeroNeighbors)
  {
    Contact contact;
    for (size_t i = 0; i < sides.size(); ++i)
      if (std::find(sides[i].otherNodes.begin(), sides[i].otherNodes.end(), zeroNeighbor) !=
              sides[i].otherNodes.end() &&
          !contains(contact.walls, wallOfSide[i]))
        contact.walls.push_back(wallOfSide[i]);
    // Along an edge between walls, the interface moves into the walls on the negative side.
    std::vector<size_t> negativeWalls;
    for (auto w : contact.walls)
      if (wallHasNegativeNeighbor[w]) negativeWalls.push_back(w);
    if (!negativeWalls.empty()) contact.walls = negativeWalls;
    contacts.push_back(contact);
  }
  if (contacts.empty()) return {};

  std::vector<size_t> commonWalls = contacts[0].walls;
  for (auto & contact : contacts)
  {
    std::vector<size_t> stillCommon;
    for (auto w : commonWalls)
      if (contains(contact.walls, w)) stillCommon.push_back(w);
    commonWalls = stillCommon;
  }
  if (!commonWalls.empty()) return surface_in_wall(commonWalls[0]);

  const size_t wall0 = contacts[0].walls[0];
  for (auto & contact : contacts)
    for (auto w : contact.walls)
      if (!contains(contacts[0].walls, w))
      {
        const stk::math::Vector3d tangent = Cross(wallAreaVector[wall0], wallAreaVector[w]);
        if (tangent.length_squared() == 0.) return {};
        return {NodalMotionConstraint::LINE, tangent.unit_vector()};
      }
  return {};
}

stk::math::Vector3d constrain_keeping_normal_speed(
    const NodalMotionConstraint & constraint, const stk::math::Vector3d & v)
{
  stk::math::Vector3d constrained;
  if (constraint.type == NodalMotionConstraint::SURFACE)
    constrained = v - Dot(v, constraint.direction) * constraint.direction;
  else if (constraint.type == NodalMotionConstraint::LINE)
    constrained = Dot(v, constraint.direction) * constraint.direction;
  else
    return v;

  // A constrained part nearly orthogonal to v would need a very large speed to keep v's normal
  // speed.
  const double vDotConstrained = Dot(v, constrained);
  const double vDotV = Dot(v, v);
  if (!(vDotConstrained > 1.e-8 * vDotV)) return v;
  return (vDotV / vDotConstrained) * constrained;
}

} // namespace krino
