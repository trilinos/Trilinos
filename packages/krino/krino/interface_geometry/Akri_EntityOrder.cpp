// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_EntityOrder.hpp>

#include <Akri_QualityMetric.hpp>

namespace krino {

bool less_by_coordinates_because_it_impacts_structure(const stk::math::Vector3d & a, const stk::math::Vector3d & b)
{
  return is_less_than_in_x_then_y_then_z(a, b);
}

bool less_by_coordinates_because_it_impacts_structure(const stk::math::Vector3d & xA,
    const stk::mesh::EntityId idA,
    const stk::math::Vector3d & xB,
    const stk::mesh::EntityId idB)
{
  if (is_less_than_in_x_then_y_then_z(xA, xB)) return true;
  if (is_less_than_in_x_then_y_then_z(xB, xA)) return false;
  return idA < idB;
}

bool less_by_global_id_because_ids_are_parallel_consistent(const stk::mesh::EntityId a, const stk::mesh::EntityId b)
{
  return a < b;
}

bool less_by_global_id_because_ids_are_parallel_consistent(
    const std::vector<stk::mesh::EntityId> & a, const std::vector<stk::mesh::EntityId> & b)
{
  return a < b;
}

bool less_conditionally_by_coordinates_because_it_impacts_structure(const bool globalIdsAreParallelConsistent,
    const stk::math::Vector3d & xA,
    const stk::mesh::EntityId idA,
    const stk::math::Vector3d & xB,
    const stk::mesh::EntityId idB)
{
  if (globalIdsAreParallelConsistent) return less_by_global_id_because_ids_are_parallel_consistent(idA, idB);
  return less_by_coordinates_because_it_impacts_structure(xA, idA, xB, idB);
}

bool less_by_global_id_because_it_doesnt_impact_structure(const stk::mesh::EntityId a, const stk::mesh::EntityId b)
{
  return a < b;
}

bool less_by_global_id_because_it_doesnt_impact_structure(
    const std::vector<stk::mesh::EntityId> & a, const std::vector<stk::mesh::EntityId> & b)
{
  return a < b;
}

bool less_by_global_id_because_it_is_requested(const stk::mesh::EntityId a, const stk::mesh::EntityId b)
{
  return a < b;
}

} // namespace krino
