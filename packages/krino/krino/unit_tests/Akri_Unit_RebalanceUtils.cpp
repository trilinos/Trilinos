// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <vector>

#include <Akri_RebalanceUtils.hpp>
#include <Akri_RebalanceUtils_Impl.hpp>
#include <Akri_Unit_Single_Element_Fixtures.hpp>
#include <Akri_Unit_MeshHelpers.hpp>

namespace krino {
namespace rebalance_utils {

namespace {

/*
 * Builds a single tri mesh with 1 level of adaptivity. Parent element is ID 1,
 * Child elem and node IDs are in the ascii art below
 *
 *            3
 *           / \
 *          / 5 \
 *       6 /_____\ 5
 *        /\     /\
 *       /  \ 3 /  \
 *      /  2 \ / 4  \
 *     /______\/_____\
 *    1       4       2
 */

void create_block_and_register_fields(SimpleStkFixture & fixture)
{
  auto & meta = fixture.meta_data();
  meta.declare_part_with_topology("block_1", stk::topology::TRIANGLE_3_2D);

  meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d> >
      (stk::topology::NODE_RANK, "coordinates");
  auto & load_field =
      meta.declare_field<stk::mesh::Field<double>>(stk::topology::ELEMENT_RANK, "element_weights");
  stk::mesh::put_field_on_mesh(load_field, meta.universal_part(), nullptr);

  fixture.commit();
}

void build_unadapted_single_tri_mesh(SimpleStkFixture & fixture)
{
  auto & bulk_data = fixture.bulk_data();
  auto & block_1 = *fixture.meta_data().get_part("block_1");

  bulk_data.modification_begin();
  if(bulk_data.parallel_rank() == 0)
  {
    stk::mesh::declare_element(bulk_data, block_1, 1, {1, 2, 3});
  }
  bulk_data.modification_end();
}

void build_one_level_adapted_single_tri_mesh(SimpleStkFixture & fixture)
{
  auto & bulk_data = fixture.bulk_data();
  auto & block_1 = *fixture.meta_data().get_part("block_1");

  bulk_data.modification_begin();
  if(bulk_data.parallel_rank() == 0)
  {
    auto parent_elem = stk::mesh::declare_element(bulk_data, block_1, 1, {1, 2, 3});
    auto child1 = stk::mesh::declare_element(bulk_data, block_1, 2, {1, 4, 6});
    auto child2 = stk::mesh::declare_element(bulk_data, block_1, 3, {4, 5, 6});
    auto child3 = stk::mesh::declare_element(bulk_data, block_1, 4, {4, 2, 5});
    auto child4 = stk::mesh::declare_element(bulk_data, block_1, 5, {6, 5, 3});
    auto family_tree = bulk_data.declare_constraint(1);
    bulk_data.declare_relation(family_tree, parent_elem, 0);
    bulk_data.declare_relation(family_tree, child1, 1);
    bulk_data.declare_relation(family_tree, child2, 2);
    bulk_data.declare_relation(family_tree, child3, 3);
    bulk_data.declare_relation(family_tree, child4, 4);
  }
  bulk_data.modification_end();
}

void build_two_level_adapted_single_tri_mesh(SimpleStkFixture & fixture)
{
  build_one_level_adapted_single_tri_mesh(fixture);

  auto & bulk_data = fixture.bulk_data();
  auto & block_1 = *fixture.meta_data().get_part("block_1");

  // Refine child elem 3 an additional time and make others into transition elements
  bulk_data.modification_begin();
  if(bulk_data.parallel_rank() == 0)
  {
    {
      auto parent = bulk_data.get_entity(stk::topology::ELEMENT_RANK, 3);
      auto child1 = stk::mesh::declare_element(bulk_data, block_1, 6, {4, 8, 7});
      auto child2 = stk::mesh::declare_element(bulk_data, block_1, 7, {7, 8, 9});
      auto child3 = stk::mesh::declare_element(bulk_data, block_1, 8, {8, 5, 9});
      auto child4 = stk::mesh::declare_element(bulk_data, block_1, 9, {9, 6, 7});

      auto family_tree = bulk_data.declare_constraint(2);
      bulk_data.declare_relation(family_tree, parent, 0);
      bulk_data.declare_relation(family_tree, child1, 1);
      bulk_data.declare_relation(family_tree, child2, 2);
      bulk_data.declare_relation(family_tree, child3, 3);
      bulk_data.declare_relation(family_tree, child4, 4);
    }

    {
      auto parent = bulk_data.get_entity(stk::topology::ELEMENT_RANK, 2);
      auto child1 = stk::mesh::declare_element(bulk_data, block_1, 10, {1, 4, 7});
      auto child2 = stk::mesh::declare_element(bulk_data, block_1, 11, {1, 7, 6});

      auto family_tree = bulk_data.declare_constraint(3);
      bulk_data.declare_relation(family_tree, parent, 0);
      bulk_data.declare_relation(family_tree, child1, 1);
      bulk_data.declare_relation(family_tree, child2, 2);
    }
    {
      auto parent = bulk_data.get_entity(stk::topology::ELEMENT_RANK, 4);
      auto child1 = stk::mesh::declare_element(bulk_data, block_1, 12, {4, 2, 8});
      auto child2 = stk::mesh::declare_element(bulk_data, block_1, 13, {2, 5, 8});

      auto family_tree = bulk_data.declare_constraint(4);
      bulk_data.declare_relation(family_tree, parent, 0);
      bulk_data.declare_relation(family_tree, child1, 1);
      bulk_data.declare_relation(family_tree, child2, 2);
    }
    {
      auto parent = bulk_data.get_entity(stk::topology::ELEMENT_RANK, 4);
      auto child1 = stk::mesh::declare_element(bulk_data, block_1, 14, {9, 5, 3});
      auto child2 = stk::mesh::declare_element(bulk_data, block_1, 15, {9, 3, 6});

      auto family_tree = bulk_data.declare_constraint(5);
      bulk_data.declare_relation(family_tree, parent, 0);
      bulk_data.declare_relation(family_tree, child1, 1);
      bulk_data.declare_relation(family_tree, child2, 2);
    }
  }
  bulk_data.modification_end();
}

}

TEST(UpdateRebalanceForAdaptivity, OneLevel)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_one_level_adapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::balance::DecompositionChangeList change_list(fixture.bulk_data(), {});
  auto parent_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);

  stk::mesh::EntityVector child_elems;
  for(int child_id=2; child_id < 6; ++child_id)
  {
    auto child_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, child_id);
    ASSERT_TRUE(mesh.is_valid(child_elem));
    child_elems.push_back(child_elem);
  }

  const int dest_proc = 2;
  change_list.set_entity_destination(parent_elem, dest_proc);

  for(auto && child_elem : child_elems)
  {
    ASSERT_FALSE(change_list.has_entity(child_elem));
  }

  impl::update_rebalance_for_adaptivity(change_list, mesh);

  for(auto && child_elem : child_elems)
  {
    EXPECT_TRUE(change_list.has_entity(child_elem));
    EXPECT_EQ(dest_proc, change_list.get_entity_destination(child_elem));
  }
}

TEST(UpdateRebalanceForAdaptivity, OneLevelChildMovedWithoutParent)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_one_level_adapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::balance::DecompositionChangeList change_list(fixture.bulk_data(), {});

  stk::mesh::EntityVector child_elems;
  for(int child_id=2; child_id < 6; ++child_id)
  {
    auto child_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, child_id);
    ASSERT_TRUE(mesh.is_valid(child_elem));
    child_elems.push_back(child_elem);
    change_list.set_entity_destination(child_elem, child_id);
  }

  impl::update_rebalance_for_adaptivity(change_list, mesh);

  for(auto && child_elem : child_elems)
  {
    EXPECT_FALSE(change_list.has_entity(child_elem));
  }
}

TEST(UpdateRebalanceForAdaptivity, TwoLevels)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_two_level_adapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::balance::DecompositionChangeList change_list(fixture.bulk_data(), {});
  auto parent_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);

  stk::mesh::EntityVector child_elems;
  for(int child_id=2; child_id < 15; ++child_id)
  {
    auto child_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, child_id);
    ASSERT_TRUE(mesh.is_valid(child_elem));
    child_elems.push_back(child_elem);
  }

  const int dest_proc = 2;
  change_list.set_entity_destination(parent_elem, dest_proc);

  impl::update_rebalance_for_adaptivity(change_list, mesh);

  for(auto && child_elem : child_elems)
  {
    EXPECT_TRUE(change_list.has_entity(child_elem));
    EXPECT_EQ(dest_proc, change_list.get_entity_destination(child_elem));
  }

  stk::mesh::EntityVector family_trees;
  mesh.get_entities(
      stk::topology::CONSTRAINT_RANK, mesh.mesh_meta_data().locally_owned_part(), family_trees);
  for(auto && ft : family_trees)
  {
    EXPECT_TRUE(change_list.has_entity(ft));
    EXPECT_EQ(dest_proc, change_list.get_entity_destination(ft));
  }
}

TEST(UpdateRebalanceForAdaptivity, TwoLevelsFirstLevelChildInitialDifferentProc)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_two_level_adapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::balance::DecompositionChangeList change_list(fixture.bulk_data(), {});
  auto parent_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);

  stk::mesh::EntityVector child_elems;
  for(int child_id=2; child_id < 15; ++child_id)
  {
    auto child_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, child_id);
    ASSERT_TRUE(mesh.is_valid(child_elem));
    child_elems.push_back(child_elem);
  }

  const int dest_proc = 2;
  change_list.set_entity_destination(parent_elem, dest_proc);
  change_list.set_entity_destination(mesh.get_entity(stk::topology::ELEMENT_RANK, 3), 5);

  impl::update_rebalance_for_adaptivity(change_list, mesh);

  for(auto && child_elem : child_elems)
  {
    EXPECT_TRUE(change_list.has_entity(child_elem));
    EXPECT_EQ(dest_proc, change_list.get_entity_destination(child_elem));
  }
}

TEST(AccumulateAdaptivityChildWeights, TwoLevelAdaptedTri)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_two_level_adapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::mesh::Field<double> & weights_field = static_cast<stk::mesh::Field<double> &>(
      *fixture.meta_data().get_field(stk::topology::ELEMENT_RANK, "element_weights"));

  auto parent_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  double & parent_weight = *stk::mesh::field_data(weights_field, parent_elem);
  parent_weight = 10.;
  stk::mesh::EntityVector child_elems;
  for (int child_id = 2; child_id < 15; ++child_id)
  {
    auto child_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, child_id);
    ASSERT_TRUE(mesh.is_valid(child_elem));
    child_elems.push_back(child_elem);

    double & weight = *stk::mesh::field_data(weights_field, child_elem);
    weight = 1.;
  }

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, weights_field);

  for (auto && child_elem : child_elems)
  {
    const double & weight = *stk::mesh::field_data(weights_field, child_elem);
    EXPECT_DOUBLE_EQ(0., weight);
  }
  EXPECT_DOUBLE_EQ(23., parent_weight);
}

TEST(AccumulateAdaptivityChildWeights, UnadaptedElement)
{
  SimpleStkFixture fixture(2, MPI_COMM_SELF);
  create_block_and_register_fields(fixture);
  build_unadapted_single_tri_mesh(fixture);

  const auto & mesh = fixture.bulk_data();

  stk::mesh::Field<double> & weights_field = static_cast<stk::mesh::Field<double> &>(
      *fixture.meta_data().get_field(stk::topology::ELEMENT_RANK, "element_weights"));

  auto parent_elem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  double & parent_weight = *stk::mesh::field_data(weights_field, parent_elem);
  parent_weight = 10.;

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, weights_field);

  EXPECT_DOUBLE_EQ(10., parent_weight);
}

TEST(Rebalance, MultipleWeightFields)
{
  SimpleStkFixture fixture(2, MPI_COMM_WORLD);

  auto & mesh = fixture.bulk_data();
  auto & meta = fixture.meta_data();

  if (mesh.parallel_size() != 2) return;

  stk::mesh::Part & block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part & block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);

  auto & coords_field = meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian2d>>(
      stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh(coords_field, meta.universal_part(), nullptr);

  auto & weights_field_1 = meta.declare_field<stk::mesh::Field<double>>(
      stk::topology::ELEMENT_RANK, "element_weights_1");
  stk::mesh::put_field_on_mesh(weights_field_1, block_1, nullptr);

  auto & weights_field_2 = meta.declare_field<stk::mesh::Field<double>>(
      stk::topology::ELEMENT_RANK, "element_weights_2");
  stk::mesh::put_field_on_mesh(weights_field_2, block_2, nullptr);

  meta.commit();

  // Create mesh with two disconnected blocks each containing two elements. Initially block_1 will
  // be
  // all owned by P0 and block_2 all by P1. Each element gets a weight of 1 in the weight field for
  // its
  // block, rebalance should put one element from each block on each proc.
  const std::vector<stk::mesh::EntityIdVector> elem_nodes{
      {1, 2, 5, 4}, {2, 3, 6, 5}, {7, 8, 11, 12}, {8, 9, 10, 11}};
  const std::vector<int> elem_procs{0, 0, 1, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_1}, {&block_2}, {&block_2}};
  std::vector<std::vector<double>> node_coords{{0, 0},
      {1, 0},
      {2, 0},
      {2, 1},
      {1, 1},
      {0, 1},
      {3, 0},
      {4, 0},
      {5, 0},
      {3, 1},
      {4, 1},
      {5, 1}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  for (auto && b : mesh.buckets(stk::topology::NODE_RANK))
  {
    for (auto && node : *b)
    {
      double * coords = stk::mesh::field_data(coords_field, node);
      coords[0] = node_coords[mesh.identifier(node) - 1][0];
      coords[1] = node_coords[mesh.identifier(node) - 1][1];
    }
  }

  stk::mesh::field_fill(1., weights_field_1);
  stk::mesh::field_fill(1., weights_field_2);

  const auto parallel_rank = mesh.parallel_rank();
  const auto & owned_part = meta.locally_owned_part();
  if (parallel_rank == 0)
  {
    EXPECT_EQ(2u,
        stk::mesh::count_selected_entities(
                  owned_part & block_1, mesh.buckets(stk::topology::ELEMENT_RANK)));
    EXPECT_EQ(0u,
        stk::mesh::count_selected_entities(
                  owned_part & block_2, mesh.buckets(stk::topology::ELEMENT_RANK)));
  }
  if (parallel_rank == 1)
  {
    EXPECT_EQ(0u,
        stk::mesh::count_selected_entities(
                  owned_part & block_1, mesh.buckets(stk::topology::ELEMENT_RANK)));
    EXPECT_EQ(2u,
        stk::mesh::count_selected_entities(
                  owned_part & block_2, mesh.buckets(stk::topology::ELEMENT_RANK)));
  }

  rebalance_utils::rebalance_mesh(
      mesh, nullptr, {weights_field_1.name(), weights_field_2.name()}, "coordinates", "parmetis");

  EXPECT_EQ(1u,
      stk::mesh::count_selected_entities(
                owned_part & block_1, mesh.buckets(stk::topology::ELEMENT_RANK)));
  EXPECT_EQ(1u,
      stk::mesh::count_selected_entities(
                owned_part & block_2, mesh.buckets(stk::topology::ELEMENT_RANK)));
}
}
}

