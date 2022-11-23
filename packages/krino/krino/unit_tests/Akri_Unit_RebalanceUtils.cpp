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
#include <Akri_RefinementInterface.hpp>
#include <Akri_AdaptivityInterface.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {
namespace rebalance_utils {

namespace {

RefinementInterface & build_refinement(stk::mesh::MetaData & meta, const bool usePercept = false)
{
  meta.enable_late_fields();

  RefinementInterface & refinement = create_refinement(meta, usePercept, sierra::Diag::sierraTimer());

  meta.disable_late_fields();

  return refinement;
}

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

  auto & coords = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_entire_mesh(coords, 2);
  auto & load_field = meta.declare_field<double>(stk::topology::ELEMENT_RANK, "element_weights");
  stk::mesh::put_field_on_mesh(load_field, meta.universal_part(), nullptr);

  fixture.commit();
}

void set_node_coordinates(const stk::mesh::BulkData & mesh, const stk::mesh::EntityId nodeId, const stk::math::Vector2d &loc)
{
  stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
  if (mesh.is_valid(node))
  {
    double* node_coords = (double*)stk::mesh::field_data(*mesh.mesh_meta_data().coordinate_field(), node);
    node_coords[0] = loc[0];
    node_coords[1] = loc[1];
  }
}


void build_unadapted_single_tri_mesh(stk::mesh::BulkData & mesh)
{
  auto & block_1 = *mesh.mesh_meta_data().get_part("block_1");
  auto & active_part = AuxMetaData::get(mesh.mesh_meta_data()).active_part();

  mesh.modification_begin();
  if(mesh.parallel_rank() == 0)
  {
    stk::mesh::PartVector elemParts{&block_1, &active_part};
    stk::mesh::declare_element(mesh, elemParts, 1, {1, 2, 3});
  }
  mesh.modification_end();
  set_node_coordinates(mesh, 1, {0.,0.});
  set_node_coordinates(mesh, 2, {1.,0.});
  set_node_coordinates(mesh, 3, {0.5,1.});
}

void build_one_level_adapted_single_tri_mesh(stk::mesh::BulkData & mesh, RefinementInterface & refinement)
{
  build_unadapted_single_tri_mesh(mesh);

  mark_elements_with_given_ids_for_refinement(mesh, refinement, {1});
  refinement.do_refinement();
}

void build_two_level_adapted_single_tri_mesh(stk::mesh::BulkData & mesh, RefinementInterface & refinement)
{
  build_one_level_adapted_single_tri_mesh(mesh, refinement);

  // Refine child elem 3 an additional time and make others into transition elements

  clear_refinement_marker(refinement);

  if(mesh.parallel_rank() == 0)
  {
    auto parentElem = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
    stk::mesh::EntityVector childElems;
    refinement.fill_children(parentElem, childElems);
    ASSERT_EQ(4u, childElems.size());

    mark_elements_for_refinement(refinement, {childElems[3]});
  }

  refinement.do_refinement();
}

}

class FixtureWithBlockAndFields : public ::testing::Test
{
public:
  FixtureWithBlockAndFields()
  {
    create_block_and_register_fields(fixture);
  }
protected:
  SimpleStkFixture fixture{2, MPI_COMM_SELF};
  stk::mesh::BulkData & mesh{fixture.bulk_data()};
};

class RebalanceForAdaptivityFixture : public FixtureWithBlockAndFields
{
public:
  RebalanceForAdaptivityFixture(bool usePercept = false) : refinement(build_refinement(fixture.meta_data(), usePercept)) {}
protected:
  stk::mesh::Entity parent_elem() const { return mesh.get_entity(stk::topology::ELEMENT_RANK, 1); }
  stk::mesh::EntityVector all_child_elems() const { stk::mesh::EntityVector allChildElems; fill_all_children(refinement, parent_elem(), allChildElems); return allChildElems; }
  void test_all_elems_have_no_destination(const stk::mesh::EntityVector & elems) const
  {
    for(auto && elem : elems)
    {
      EXPECT_FALSE(change_list.has_entity(elem));
    }
  }
  void test_all_entities_have_destination(const stk::mesh::EntityVector & entities, const int dest_proc) const
  {
    for(auto && entity : entities)
    {
      EXPECT_TRUE(change_list.has_entity(entity));
      EXPECT_EQ(dest_proc, change_list.get_entity_destination(entity));
    }
  }
  void test_all_family_tree_entities_have_destination(const int dest_proc) const
  {
    stk::mesh::EntityVector family_trees;
    mesh.get_entities(
        stk::topology::CONSTRAINT_RANK, mesh.mesh_meta_data().locally_owned_part(), family_trees);
    for(auto && ft : family_trees)
    {
      if (mesh.num_elements(ft) > 0)
      {
        EXPECT_TRUE(change_list.has_entity(ft));
        EXPECT_EQ(dest_proc, change_list.get_entity_destination(ft));
      }
    }
  }
  stk::mesh::Field<double> & weights_field() const
  {
    stk::mesh::Field<double> & weightsField = static_cast<stk::mesh::Field<double> &>(
      *fixture.meta_data().get_field(stk::topology::ELEMENT_RANK, "element_weights"));
    return weightsField;
  }
  void set_weight_for_elem(stk::mesh::Entity elem, const double newWeight) const
  {
    stk::mesh::Field<double> & weightsField = weights_field();
    double & weight = *stk::mesh::field_data(weightsField, elem);
    weight = newWeight;
  }
  void test_element_has_weight(const stk::mesh::Entity & elem, const double goldWeight) const
  {
    stk::mesh::Field<double> & weightsField = weights_field();
    const double & weight = *stk::mesh::field_data(weightsField, elem);
    EXPECT_DOUBLE_EQ(goldWeight, weight);
  }
  void test_elements_have_weight(const stk::mesh::EntityVector & elems, const double goldWeight) const
  {
    stk::mesh::Field<double> & weightsField = weights_field();
    for (auto && elem : elems)
    {
      const double & weight = *stk::mesh::field_data(weightsField, elem);
      EXPECT_DOUBLE_EQ(goldWeight, weight);
    }
  }
  RefinementInterface & refinement;
  stk::balance::DecompositionChangeList change_list{fixture.bulk_data(), {}};
};

class RebalanceForAdaptivityFixtureForPercept : public RebalanceForAdaptivityFixture
{
public:
  RebalanceForAdaptivityFixtureForPercept() : RebalanceForAdaptivityFixture(true) {}
private:
};

TEST_F(RebalanceForAdaptivityFixture, OneLevel)
{
  build_one_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector childElems = all_child_elems();

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);

  test_all_elems_have_no_destination(childElems);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(childElems, destProc);
}

TEST_F(RebalanceForAdaptivityFixture, OneLevelChildMovedWithoutParent)
{
  build_one_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector childElems = all_child_elems();

  test_all_elems_have_no_destination(childElems);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_elems_have_no_destination(childElems);
}

TEST_F(RebalanceForAdaptivityFixture, TwoLevels)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(allChildElems, destProc);
  test_all_family_tree_entities_have_destination(destProc);
}

TEST_F(RebalanceForAdaptivityFixture, TwoLevelsFirstLevelChildInitialDifferentProc)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);
  change_list.set_entity_destination(allChildElems[2], 5);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(allChildElems, destProc);
  test_all_family_tree_entities_have_destination(destProc);
}

TEST_F(RebalanceForAdaptivityFixture, AccumulateAdaptivityChildWeightsTwoLevelAdaptedTri)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());
  const stk::mesh::EntityVector allButLastChildElems(allChildElems.begin(), allChildElems.begin()+13);

  for (auto childElem : allButLastChildElems)
    set_weight_for_elem(childElem, 1.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  test_elements_have_weight(allButLastChildElems, 0.);
  test_element_has_weight(parent_elem(), 23.);
}

TEST_F(RebalanceForAdaptivityFixture, AccumulateAdaptivityChildWeightsUnadaptedElement)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  test_element_has_weight(parent_elem(), 10.);
}

// FIXME: Temporary version of tests for percept
#if 1
TEST_F(RebalanceForAdaptivityFixtureForPercept, OneLevel)
{
  build_one_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector childElems = all_child_elems();

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);

  test_all_elems_have_no_destination(childElems);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(childElems, destProc);
}

TEST_F(RebalanceForAdaptivityFixtureForPercept, OneLevelChildMovedWithoutParent)
{
  build_one_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector childElems = all_child_elems();

  test_all_elems_have_no_destination(childElems);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_elems_have_no_destination(childElems);
}

TEST_F(RebalanceForAdaptivityFixtureForPercept, TwoLevels)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(allChildElems, destProc);
  test_all_family_tree_entities_have_destination(destProc);
}

TEST_F(RebalanceForAdaptivityFixtureForPercept, TwoLevelsFirstLevelChildInitialDifferentProc)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  const int destProc = 2;
  change_list.set_entity_destination(parent_elem(), destProc);
  change_list.set_entity_destination(allChildElems[2], 5);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  test_all_entities_have_destination(allChildElems, destProc);
  test_all_family_tree_entities_have_destination(destProc);
}

TEST_F(RebalanceForAdaptivityFixtureForPercept, AccumulateAdaptivityChildWeightsTwoLevelAdaptedTri)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());
  const stk::mesh::EntityVector allButLastChildElems(allChildElems.begin(), allChildElems.begin()+13);

  for (auto childElem : allButLastChildElems)
    set_weight_for_elem(childElem, 1.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  test_elements_have_weight(allButLastChildElems, 0.);
  test_element_has_weight(parent_elem(), 23.);
}

TEST_F(RebalanceForAdaptivityFixtureForPercept, AccumulateAdaptivityChildWeightsUnadaptedElement)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  test_element_has_weight(parent_elem(), 10.);
}
#endif

TEST(Rebalance, MultipleWeightFields)
{
  SimpleStkFixture fixture(2, MPI_COMM_WORLD);

  auto & mesh = fixture.bulk_data();
  auto & meta = fixture.meta_data();

  if (mesh.parallel_size() != 2) return;

  stk::mesh::Part & block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part & block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);

  auto & coords_field = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh(coords_field, meta.universal_part(), 2, nullptr);

  auto & weights_field_1 = meta.declare_field<double>(stk::topology::ELEMENT_RANK, "element_weights_1");
  stk::mesh::put_field_on_mesh(weights_field_1, block_1, nullptr);

  auto & weights_field_2 = meta.declare_field<double>(stk::topology::ELEMENT_RANK, "element_weights_2");
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

  rebalance_utils::rebalance_mesh(mesh,
      nullptr,
      nullptr,
      {weights_field_1.name(), weights_field_2.name()},
      "coordinates",
      10,
      "parmetis");

  EXPECT_EQ(1u,
      stk::mesh::count_selected_entities(
                owned_part & block_1, mesh.buckets(stk::topology::ELEMENT_RANK)));
  EXPECT_EQ(1u,
      stk::mesh::count_selected_entities(
                owned_part & block_2, mesh.buckets(stk::topology::ELEMENT_RANK)));
}
}
}

