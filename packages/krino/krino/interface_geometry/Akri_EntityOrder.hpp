// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_ENTITYORDER_HPP_
#define KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_ENTITYORDER_HPP_

#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

// The orders in which krino compares mesh nodes, points, and (Akri_SubElementNodeOrder.hpp)
// sub-element nodes, each named for why it is used:
// - less_by_coordinates_because_it_impacts_structure: for choices that change a decomposition (its
//   nodes, elements or phases).  By location in single precision, x then y then z, so the same for
//   any ids and any parallel decomposition; by id only where locations tie exactly.
// - less_conditionally_by_coordinates_because_it_impacts_structure: for such choices made by id
//   (less_by_global_id_because_ids_are_parallel_consistent) where the global ids are parallel
//   consistent (CDFEM_Support::get_global_ids_are_parallel_consistent: the same mesh has the same ids
//   however it is decomposed), and by location otherwise (refinement, which numbers new entities
//   differently for different decompositions and builds, clears it).  Where the keys of the two
//   orders are costly, the caller picks between them on the flag.
// - less_by_global_id_because_it_doesnt_impact_structure: for choices that change only numbering or
//   storage.
// - less_by_global_id_because_it_is_requested: for the methods that ask for ids
//   (CUT_QUADS_BY_GLOBAL_IDENTIFIER), although they change the decomposition.
// Every comparison of these uses one of them.

bool less_by_coordinates_because_it_impacts_structure(const stk::math::Vector3d & a, const stk::math::Vector3d & b);
bool less_by_coordinates_because_it_impacts_structure(const stk::math::Vector3d & xA,
    const stk::mesh::EntityId idA,
    const stk::math::Vector3d & xB,
    const stk::mesh::EntityId idB);

// For std::sort and the like, which cannot take the overloaded name.
struct LessByCoordinatesBecauseItImpactsStructure
{
  bool operator()(const stk::math::Vector3d & a, const stk::math::Vector3d & b) const
  {
    return less_by_coordinates_because_it_impacts_structure(a, b);
  }
};

bool less_by_global_id_because_ids_are_parallel_consistent(const stk::mesh::EntityId a, const stk::mesh::EntityId b);
// Lexicographically.
bool less_by_global_id_because_ids_are_parallel_consistent(
    const std::vector<stk::mesh::EntityId> & a, const std::vector<stk::mesh::EntityId> & b);

bool less_conditionally_by_coordinates_because_it_impacts_structure(const bool globalIdsAreParallelConsistent,
    const stk::math::Vector3d & xA,
    const stk::mesh::EntityId idA,
    const stk::math::Vector3d & xB,
    const stk::mesh::EntityId idB);

bool less_by_global_id_because_it_doesnt_impact_structure(const stk::mesh::EntityId a, const stk::mesh::EntityId b);
// Lexicographically.
bool less_by_global_id_because_it_doesnt_impact_structure(
    const std::vector<stk::mesh::EntityId> & a, const std::vector<stk::mesh::EntityId> & b);

bool less_by_global_id_because_it_is_requested(const stk::mesh::EntityId a, const stk::mesh::EntityId b);

} // namespace krino

#endif /* KRINO_KRINO_INTERFACE_GEOMETRY_AKRI_ENTITYORDER_HPP_ */
