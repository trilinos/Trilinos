// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <vector>

#include <gtest/gtest.h>

#include <Akri_EntityOrder.hpp>

namespace krino
{

TEST(EntityOrder, byCoordinates_xThenYThenZInSinglePrecisionThenId)
{
  const stk::math::Vector3d a(0., 1., 0.);
  const stk::math::Vector3d b(1., 0., 0.);
  EXPECT_TRUE(less_by_coordinates_because_it_impacts_structure(a, b));
  EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(b, a));
  // Equal in single precision: by id.
  const stk::math::Vector3d c(0.5, 1., 0.);
  const stk::math::Vector3d cNudged(0.5 + 1.e-12, 1., 0.);
  EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(c, cNudged));
  EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(cNudged, c));
  EXPECT_TRUE(less_by_coordinates_because_it_impacts_structure(cNudged, 1, c, 2));
  EXPECT_FALSE(less_by_coordinates_because_it_impacts_structure(cNudged, 2, c, 1));
  EXPECT_TRUE(less_by_coordinates_because_it_impacts_structure(a, 2, b, 1));
}

TEST(EntityOrder, conditionallyByCoordinates_byIdWhereIdsAreParallelConsistent)
{
  const stk::math::Vector3d a(0., 1., 0.);
  const stk::math::Vector3d b(1., 0., 0.);
  EXPECT_TRUE(less_conditionally_by_coordinates_because_it_impacts_structure(true, b, 1, a, 2));
  EXPECT_FALSE(less_conditionally_by_coordinates_because_it_impacts_structure(false, b, 1, a, 2));
  EXPECT_TRUE(less_conditionally_by_coordinates_because_it_impacts_structure(false, a, 2, b, 1));

}

TEST(EntityOrder, byGlobalIdBecauseIdsAreParallelConsistent_lexicographicForLists)
{
  EXPECT_TRUE(less_by_global_id_because_ids_are_parallel_consistent(stk::mesh::EntityId(1), stk::mesh::EntityId(2)));
  EXPECT_TRUE(less_by_global_id_because_ids_are_parallel_consistent(
      std::vector<stk::mesh::EntityId>{1, 2}, std::vector<stk::mesh::EntityId>{1, 3}));
  EXPECT_FALSE(less_by_global_id_because_ids_are_parallel_consistent(
      std::vector<stk::mesh::EntityId>{1, 3}, std::vector<stk::mesh::EntityId>{1, 2}));
}

TEST(EntityOrder, byGlobalId)
{
  EXPECT_TRUE(less_by_global_id_because_it_doesnt_impact_structure(stk::mesh::EntityId(1), stk::mesh::EntityId(2)));
  EXPECT_FALSE(less_by_global_id_because_it_doesnt_impact_structure(stk::mesh::EntityId(2), stk::mesh::EntityId(2)));
  EXPECT_TRUE(less_by_global_id_because_it_doesnt_impact_structure(
      std::vector<stk::mesh::EntityId>{1, 2}, std::vector<stk::mesh::EntityId>{1, 3}));
  EXPECT_TRUE(less_by_global_id_because_it_is_requested(stk::mesh::EntityId(1), stk::mesh::EntityId(2)));
}

} // namespace krino
