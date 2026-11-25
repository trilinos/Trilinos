// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <memory>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

#include <Akri_MeshHelpers.hpp>
#include <Akri_Unit_MeshHelpers.hpp>

namespace krino {

namespace {
void test_and_cleanup_internal_side(stk::mesh::BulkData & mesh, const stk::mesh::Part & block_1, const stk::mesh::Part & block_2)
{
  EXPECT_TRUE(check_face_and_edge_relations(mesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_TRUE(check_shared_entity_nodes(mesh));

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_entities( mesh, mesh.mesh_meta_data().side_rank(), sides );

  EXPECT_TRUE(1 == sides.size());

  for (auto && side : sides)
  {
    EXPECT_TRUE(mesh.bucket(side).member(block_1));
    EXPECT_TRUE(mesh.bucket(side).member(block_2));
  }

  // cleanup
  mesh.modification_begin();
  for (auto && side : sides)
  {
    EXPECT_TRUE(disconnect_and_destroy_entity(mesh, side));
  }
  mesh.modification_end();
}
}

auto create_2D_mesh(const stk::ParallelMachine & pm)
{
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(pm).set_spatial_dimension(2).create();
  return bulk;
}

TEST(MeshHelpers, DeclareElementSide)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  /*  Mesh
   *    4---5---6  P0 owns nodes 1,2,4,5; P, elem 1
   *    | 1 | 2 |  P1 : 3,6, elem 2
   *    1---2---3
   */


  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& surface_1 = meta.declare_part_with_topology("surface_1", stk::topology::LINE_2);

  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 5, 4}, {2, 3, 6, 5} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_2}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  // All elements and nodes appear everywhere via aura
  stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
  stk::mesh::Entity node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  stk::mesh::Entity node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

  std::vector<stk::mesh::Entity> side_nodes;
  side_nodes.push_back(node2);
  side_nodes.push_back(node5);

  mesh.initialize_face_adjacent_element_graph();

  {
    // Case 1: Declare side on owning processor for both elements
    mesh.modification_begin();

    if (mesh.bucket(element1).owned())
    {
      const unsigned elem1_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element1, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element1, elem1_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }
    if (mesh.bucket(element2).owned())
    {
      const unsigned elem2_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element2, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element2, elem2_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }

  {
    // Case 2: Declare side of element1 on processor that owns element1
    mesh.modification_begin();

    if (mesh.bucket(element1).owned())
    {
      const unsigned elem1_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element1, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element1, elem1_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }

  {
    // Case 3: Declare side of element2 on processor that owns element2
    mesh.modification_begin();

    if (mesh.bucket(element2).owned())
    {
      const unsigned elem2_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element2, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element2, elem2_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }
}

TEST(MeshHelpers, FullyCoincidentVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // This test will create a two element mesh (quad4 elements), on more than 2 processors only
  // ranks 0 and 1 will have any elements. We test larger number of processors to ensure that
  // we get a parallel-consistent result to avoid potential parallel hangs in the full app.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 3, 4}, {2, 3, 4, 1} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_1}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_FALSE(ok);
}

TEST(MeshHelpers, PartiallyCoincidentActiveVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 3, 4}, {1, 2, 5, 6} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1, &active_part}, {&block_1, &active_part}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_FALSE(ok);
}

TEST(MeshHelpers, NotCoincidentActiveDegenerateVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 2, 3}, {3, 2, 2, 4} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1, &active_part}, {&block_1, &active_part}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_TRUE(ok);
}

TEST(MeshHelpers, EdgeAttachedToOnlyOneOfTheElementsThatHaveIt_failsFullConnectivityCheckOnly)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(pm) != 1) return;

  // Two tets share the face 1-2-3.  The edge 1-2 is an edge of both, but is attached only to
  // element 1.
  std::unique_ptr<stk::mesh::BulkData> meshPtr =
      stk::mesh::MeshBuilder(pm).set_spatial_dimension(3).create();
  stk::mesh::BulkData & mesh = *meshPtr;
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::Part & block_1 =
      meta.declare_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  stk::mesh::Part & edge_part = meta.declare_part_with_topology("edge_part", stk::topology::LINE_2);
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{{1, 2, 3, 4}, {1, 3, 2, 5}};
  const std::vector<int> elem_procs{0, 0};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_1}};
  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  mesh.modification_begin();
  stk::mesh::Entity edge = mesh.declare_edge(1, stk::mesh::PartVector{&edge_part});
  mesh.declare_relation(edge, mesh.get_entity(stk::topology::NODE_RANK, 1), 0);
  mesh.declare_relation(edge, mesh.get_entity(stk::topology::NODE_RANK, 2), 1);
  attach_entity_to_element(mesh, edge, element1);
  mesh.modification_end();

  EXPECT_FALSE(check_face_and_edge_relations(mesh, true));
  EXPECT_TRUE(check_face_and_edge_relations(mesh, false));

  // Attached to element 1 with the wrong ordinal: incorrect even if not fully connected.
  const auto [edgeOrdinal, edgePermutation] =
      determine_ordinal_and_permutation(mesh, element1, edge);
  const stk::mesh::ConnectivityOrdinal wrongEdgeOrdinal = (edgeOrdinal + 1) % 6;
  mesh.modification_begin();
  mesh.destroy_relation(element1, edge, edgeOrdinal);
  mesh.declare_relation(element1, edge, wrongEdgeOrdinal, edgePermutation);
  mesh.modification_end();

  EXPECT_FALSE(check_face_and_edge_relations(mesh, false));

  // Not attached to any element.
  mesh.modification_begin();
  mesh.destroy_relation(element1, edge, wrongEdgeOrdinal);
  mesh.modification_end();

  EXPECT_FALSE(check_face_and_edge_relations(mesh, false));
}

TEST(MeshHelpers,
    SharedSideOwnedByProcWithoutItsElements_withoutAura_ownershipFixedToProcWithElement)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(pm) != 2) return;

  /*  Mesh: elements 1, 3, 4 on P0 and element 2 on P1.  Elements 1 and 2 share the side 2-5, owned
   *  by P0.  Then element 1 moves to P1, so P0 owns the side but none of its elements, while
   *  elements 3 and 4 keep its nodes on P0.
   *    9--10
   *    | 4 |
   *    4---5---6
   *    | 1 | 2 |
   *    1---2---3
   *    | 3 |
   *    7---8
   */
  std::unique_ptr<stk::mesh::BulkData> meshPtr =
      stk::mesh::MeshBuilder(pm)
          .set_spatial_dimension(2)
          .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
          .create();
  stk::mesh::BulkData & mesh = *meshPtr;
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Part & block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part & surface_1 = meta.declare_part_with_topology("surface_1", stk::topology::LINE_2);
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{
      {1, 2, 5, 4}, {2, 3, 6, 5}, {7, 8, 2, 1}, {4, 5, 10, 9}};
  const std::vector<int> elem_procs{0, 1, 0, 0};
  std::vector<stk::mesh::PartVector> elem_parts(4, {&block_1});
  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const std::vector<stk::mesh::Entity> sideNodes{
      mesh.get_entity(stk::topology::NODE_RANK, 2), mesh.get_entity(stk::topology::NODE_RANK, 5)};
  mesh.initialize_face_adjacent_element_graph();
  mesh.modification_begin();
  for (stk::mesh::EntityId elemId : {1, 2})
  {
    const stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
    if (mesh.is_valid(elem) && mesh.bucket(elem).owned())
    {
      const unsigned sideOrdinal = stk::mesh::get_entity_subcell_id(
          mesh, elem, meta.side_rank(), surface_1.topology(), sideNodes);
      mesh.declare_element_side(elem, sideOrdinal, stk::mesh::ConstPartVector{&surface_1});
    }
  }
  mesh.modification_end();

  // Whether the one side is owned by the given processor wherever it is (a processor without any of
  // its elements can drop it).  Collective, so that a failure returns from the test on every
  // processor.
  auto is_side_owned_by_on_all_procs = [&](const int owner)
  {
    std::vector<stk::mesh::Entity> sides;
    stk::mesh::get_entities(mesh, meta.side_rank(), sides);
    const bool isOwnedByOwnerHere = sides.size() <= 1 &&
        (sides.empty() || mesh.parallel_owner_rank(sides[0]) == owner) &&
        (mesh.parallel_rank() != owner || sides.size() == 1);
    return stk::is_true_on_all_procs(pm, isOwnedByOwnerHere);
  };

  ASSERT_TRUE(is_side_owned_by_on_all_procs(0));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));

  std::vector<stk::mesh::EntityProc> elemsToMove;
  const stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  if (mesh.is_valid(element1) && mesh.bucket(element1).owned())
    elemsToMove.emplace_back(element1, 1);
  mesh.change_entity_owner(elemsToMove);

  ASSERT_TRUE(is_side_owned_by_on_all_procs(0));
  EXPECT_FALSE(check_face_and_edge_ownership(mesh));

  EXPECT_TRUE(fix_face_and_edge_ownership(mesh));

  ASSERT_TRUE(is_side_owned_by_on_all_procs(1));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_FALSE(fix_face_and_edge_ownership(mesh));
}
}
