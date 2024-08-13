// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_io/FillMesh.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <gtest/gtest.h>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
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
#include <stk_util/diag/Timer.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {
namespace rebalance_utils {

namespace {

RefinementInterface & build_refinement(stk::mesh::MetaData & meta)
{
  meta.enable_late_fields();

  RefinementInterface & refinement = KrinoRefinement::create(meta);

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
  FixtureWithBlockAndFields(MPI_Comm comm) :
    fixture(2, comm),
    mesh(fixture.bulk_data())
  {
    create_block_and_register_fields(fixture);
  }
protected:
  SimpleStkFixture fixture;
  stk::mesh::BulkData & mesh;
};

class RebalanceForAdaptivityFixture : public FixtureWithBlockAndFields
{
public:
  RebalanceForAdaptivityFixture(MPI_Comm comm = MPI_COMM_SELF) :
    FixtureWithBlockAndFields(comm),
    refinement(build_refinement(fixture.meta_data())) {}
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

class ParallelRebalanceForAdaptivityFixture3D : public ::testing::Test
{
public:
  ParallelRebalanceForAdaptivityFixture3D() :
    bulk(stk::mesh::MeshBuilder(MPI_COMM_WORLD)
      .set_spatial_dimension(3)
      .create()),
    meta(bulk->mesh_meta_data_ptr()),
    load_field(meta->declare_field<double>(stk::topology::ELEMENT_RANK, "element_weights")),
    change_list(*bulk, {})
  {
    AuxMetaData::create(*meta);
    double zero_val = 0.;
    stk::mesh::put_field_on_mesh(load_field, meta->universal_part(), &zero_val);

    stk::io::fill_mesh("generated:2x2x4|bbox:-1,-1,-1,1,1,1|tets", *bulk);
    refinement_ptr = &build_refinement(*meta);
  }

  void refine_z_boundary()
  {
    RefinementInterface & refinement = *refinement_ptr;
    clear_refinement_marker(refinement);
    auto elem_buckets = 
      bulk->get_buckets(stk::topology::ELEMENT_RANK, bulk->mesh_meta_data().locally_owned_part());
    for(auto && bucket : elem_buckets)
    {
      for(auto && elem : *bucket)
      {
        if(refinement.is_parent(elem)) continue;
        stk::mesh::Entity const * entity_node_rels = bulk->begin_nodes(elem);
        int num_entity_nodes = bulk->num_nodes(elem);
        for(int ni = 0; ni < num_entity_nodes; ++ni)
        {
          stk::mesh::Entity node = entity_node_rels[ni];
          auto * coords = field_data<double>(*meta->coordinate_field(), node);
          if(std::fabs(coords[meta->spatial_dimension()-1]-1.) <= 1e-6)
          {
            int * elemMarker = field_data<int>(refinement.get_marker_field_and_sync_to_host(), elem);
            *elemMarker = static_cast<int>(Refinement::RefinementMarker::REFINE);
          }
        }
      }
    }
    refinement.do_refinement();
  }

  void set_active_element_weights()
  {
    auto & refinement = *refinement_ptr;
    auto elem_buckets = 
      bulk->get_buckets(stk::topology::ELEMENT_RANK, bulk->mesh_meta_data().locally_owned_part());
    for(auto && bucket : elem_buckets)
    {
      for(auto && elem : *bucket)
      {
        auto * elem_load = field_data<double>(load_field, elem);
        if(refinement.is_parent(elem)) *elem_load = 0.;
        else *elem_load = 1.;
      }
    }
  }

  void execute_rebalance()
  {
    rebalance_mesh(*bulk, refinement_ptr, nullptr, 
      load_field.name(), 
      meta->coordinate_field()->name(), 
      {meta->universal_part()}, 10, "rcb");
  }
protected:
  std::unique_ptr<stk::mesh::BulkData> bulk;
  std::shared_ptr<stk::mesh::MetaData> meta;
  stk::mesh::Field<double> & load_field;
  RefinementInterface * refinement_ptr;
  stk::balance::DecompositionChangeList change_list;
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

TEST_F(RebalanceForAdaptivityFixture, TwoLevelsLeafChildInitialDifferentProc)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  //Leaf children must move with their immediate parents. Build set of parents of leaf children and assign
  //each a proc to move to
  stk::mesh::EntityVector leaf_children;
  std::set<stk::mesh::Entity> leaf_parents;
  fill_leaf_children(refinement, parent_elem(), leaf_children);
  for (auto && leaf_child : leaf_children)
  {
    leaf_parents.insert(refinement.get_parent(leaf_child));
    change_list.set_entity_destination(leaf_child, 10);
  }
  int destProc = 2;
  for (auto && leaf_parent : leaf_parents)
  {
    change_list.set_entity_destination(leaf_parent, destProc);
    destProc++;
  }
  //The grandparent should not be allowed to move if any of its children are being moved. Give it
  //a destination to ensure it gets zeroed out
  change_list.set_entity_destination(parent_elem(), destProc);

  impl::update_rebalance_for_adaptivity(change_list, refinement, mesh);

  //Verify all leaf children are moving with their parents
  for (auto && leaf_child : leaf_children)
  {
    auto parent = refinement.get_parent(leaf_child);
    EXPECT_TRUE(change_list.has_entity(leaf_child));
    EXPECT_TRUE(change_list.has_entity(parent));
    EXPECT_EQ(change_list.get_entity_destination(parent), change_list.get_entity_destination(leaf_child));
  }
  //Verify grandparent is staying put
  EXPECT_FALSE(change_list.has_entity(parent_elem()));
}

TEST_F(RebalanceForAdaptivityFixture, AccumulateAdaptivityChildWeightsTwoLevelAdaptedTri)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  const stk::mesh::EntityVector allChildElems = all_child_elems();
  ASSERT_EQ(14u, allChildElems.size());

  for (auto childElem : allChildElems)
    set_weight_for_elem(childElem, 1.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  //Only parents of leaf children should have weights accumulated. All other weights untouched
  //Build set of parents of leaf children
  stk::mesh::EntityVector leaf_children;
  std::set<stk::mesh::Entity> leaf_parents;
  fill_leaf_children(refinement, parent_elem(), leaf_children);
  for (auto && leaf_child : leaf_children)
  {
    leaf_parents.insert(refinement.get_parent(leaf_child));
  }

  //Test that all elements not in the set of parents of leaf children have original weight and those in the set
  //have weight equal to their number of children + their weight
  test_element_has_weight(parent_elem(), 10.);
  for (auto && childElem : allChildElems )
  {
    if(leaf_parents.find(childElem) != leaf_parents.end())
      test_element_has_weight(childElem, refinement.get_num_children(childElem) + 1);
    else if(!refinement.is_parent(childElem))
      test_element_has_weight(childElem, 0.);
    else
      test_element_has_weight(childElem, 1.);
  }
}

TEST_F(RebalanceForAdaptivityFixture, AccumulateAdaptivityChildWeightsUnadaptedElement)
{
  build_two_level_adapted_single_tri_mesh(mesh, refinement);

  set_weight_for_elem(parent_elem(), 10.);

  impl::accumulate_adaptivity_child_weights_to_parents(mesh, refinement, weights_field());

  test_element_has_weight(parent_elem(), 10.);
}

TEST_F(ParallelRebalanceForAdaptivityFixture3D, TwoLevelsLeafChildenConstrainedUnadaptedFree)
{
  auto & refinement = *refinement_ptr;
  for(unsigned i=0; i<2; i++)
  {
    refine_z_boundary();
  }

  //Leaf children must move with their immediate parents. Build set of parents of leaf children and assign
  //each a proc to move to
  std::set<stk::mesh::Entity> free_parents;
  stk::mesh::EntityVector children;
  auto elem_buckets = 
    bulk->get_buckets(stk::topology::ELEMENT_RANK, meta->locally_owned_part());
  int free_destProc = 10000;
  int destProc = 2;
  for(auto && bucket : elem_buckets)
  {
    for(auto && elem : *bucket)
    {
      if(refinement.is_parent(elem))
      {
        children.clear();
        refinement.fill_children(elem, children);
        if(std::any_of(children.cbegin(), children.cend(), 
          [&](stk::mesh::Entity e){ return refinement.is_parent(e); }))
        {
          //expectation is that these moves will be deleted
          change_list.set_entity_destination(elem, free_destProc);
        }
        else 
        {
          //These moves should trigger children to move along with them
          change_list.set_entity_destination(elem, destProc);
          free_parents.insert(elem);
          destProc++;
        }
      }
      //These moves should be allowed for unadapted elements, but not allowed for leaf elements
      else change_list.set_entity_destination(elem, free_destProc);
    }
  }

  impl::update_rebalance_for_adaptivity(change_list, refinement, *bulk);

  //Verify all leaf children are moving with their parents and unadapted elements are free to move
  for(auto && bucket : elem_buckets)
  {
    for(auto && elem : *bucket)
    {
      if(refinement.is_parent(elem))
      {
        if(free_parents.find(elem) == free_parents.end())
        {
          //Parent has children that are also parent, isn't allowed to move
          EXPECT_FALSE(change_list.has_entity(elem));
        }
      }
      else if(refinement.is_child(elem))
      {
        auto parent = refinement.get_parent(elem);
        if(free_parents.find(parent) != free_parents.end())
        {
          //Leaf element has parent whose children are all leaf elements, they can move, but need to do so together
          EXPECT_TRUE(change_list.has_entity(elem));
          EXPECT_TRUE(change_list.has_entity(parent));
          EXPECT_EQ(change_list.get_entity_destination(parent), change_list.get_entity_destination(elem));
        }
        //otherwise, leaf element has a parent with children who are also parent. Should not move
        else EXPECT_FALSE(change_list.has_entity(elem));
      }
      //These are unadapted elements, can move freely
      else EXPECT_EQ(change_list.get_entity_destination(elem), free_destProc);
    }
  }
}

TEST_F(ParallelRebalanceForAdaptivityFixture3D, ParentChildRebalanceRules)
{
  auto & refinement = *refinement_ptr;
  for(unsigned i=0; i<3; i++)
  {
    refine_z_boundary();
  }

  set_active_element_weights();

  execute_rebalance();

  //Verify leaf elements have their parent element on the same processor and flag for coarsening
  {
    clear_refinement_marker(refinement);
    FieldRef markerField = refinement.get_marker_field_and_sync_to_host();
    auto elem_buckets = 
      bulk->get_buckets(stk::topology::ELEMENT_RANK, meta->locally_owned_part());
    for(auto && bucket : elem_buckets)
    {
      for(auto && elem : *bucket)
      {
        if(refinement.is_parent(elem) || (!refinement.is_child(elem))) continue;
        auto parent = refinement.get_parent(elem);
        EXPECT_TRUE(parent != stk::mesh::Entity());
        EXPECT_TRUE(refinement.is_parent(parent));
        int * elemMarker = field_data<int>(markerField, elem);
        *elemMarker = static_cast<int>(Refinement::RefinementMarker::COARSEN);
      }
    }
    refinement.do_refinement();
  }

  //Verify leaf elements still have their parent element on the same processor after coarsening
  {
    auto elem_buckets = 
      bulk->get_buckets(stk::topology::ELEMENT_RANK, meta->locally_owned_part());
    for(auto && bucket : elem_buckets)
    {
      for(auto && elem : *bucket)
      {
        if(refinement.is_parent(elem) || (!refinement.is_child(elem))) continue;
        auto parent = refinement.get_parent(elem);
        EXPECT_TRUE(parent != stk::mesh::Entity());
        EXPECT_TRUE(refinement.is_parent(parent));
      }
    }
  }
}

TEST_F(ParallelRebalanceForAdaptivityFixture3D, EffectiveLoadBalance)
{
  auto & refinement = *refinement_ptr;

  for(unsigned i=0; i<5; i++)
  {
    refine_z_boundary();
    set_active_element_weights();
    execute_rebalance();
  }

  unsigned elem_cnt = 0;
  auto elem_buckets = 
    bulk->get_buckets(stk::topology::ELEMENT_RANK, bulk->mesh_meta_data().locally_owned_part());
  for(auto && bucket : elem_buckets)
  {
    for(auto && elem : *bucket)
    {
      if(!refinement.is_parent(elem)) elem_cnt++;
    }
  }
  
  double elem_max = stk::get_global_max(bulk->parallel(), elem_cnt);
  double elem_min = stk::get_global_min(bulk->parallel(), elem_cnt);

  EXPECT_LT(elem_max/elem_min, 1.2);
}

TEST(Rebalance, MultipleWeightFields)
{
  if (!rebalance_utils::have_parmetis())
    return;

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

